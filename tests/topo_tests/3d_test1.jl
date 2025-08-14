using StaticArrays, WriteVTK
using OrderedCollections, Bumper
using LinearAlgebra, Statistics
using SmallCollections, Chairmarks

# Include the main module (professional approach)
include("../../src/small_utils.jl")
include("../../src/lagrange_utils.jl")
include("../../src/topo.jl")
include("../../src/vtkexports.jl")
include("../../src/element_refinement.jl")
include("../../src/element_coarsening.jl")
include("../../src/triangulation.jl")
# =============================================================================
# 3D Topology Test - Experimental Code
# =============================================================================

# Grid parameters
# nx = 30; ny = 30; nz = 30
nx,ny,nz = 10,10,4
dx = 1.0; dy = 1.0; dz = 1.0

# Generate coordinate ranges
x_coords = range(0, nx*dx, length=nx+1)
y_coords = range(0, ny*dy, length=ny+1)
z_coords = range(0, nz*dz, length=nz+1)

# Create coordinate points
coords = [SA[x,y,z] for x in x_coords, y in y_coords, z in z_coords]

# Initialize topology
topo = Topology{3}()

# Add all nodes to topology
add_node!.(coords, Ref(topo))

# Create linear indices for easy access
idxs = LinearIndices((nx+1, ny+1, nz+1))

# Test coordinate access
println("Coordinate at (2,2,2): ", coords[idxs[2,2,2]])

# =============================================================================
# Create hexahedral elements
# =============================================================================

for I in CartesianIndices((nx, ny, nz))
    i, j, k = Tuple(I)
    
    # Define face IDs for each hexahedron
    _face_ids = [
        [idxs[i,j,k], idxs[i+1,j,k], idxs[i+1,j+1,k], idxs[i,j+1,k]],           # bottom
        [idxs[i,j,k+1], idxs[i+1,j,k+1], idxs[i+1,j+1,k+1], idxs[i,j+1,k+1]],   # top
        [idxs[i,j,k], idxs[i+1,j,k], idxs[i+1,j,k+1], idxs[i,j,k+1]],           # front 
        [idxs[i+1,j,k], idxs[i+1,j+1,k], idxs[i+1,j+1,k+1], idxs[i+1,j,k+1]],   # right 
        [idxs[i,j,k], idxs[i,j,k+1], idxs[i,j+1,k+1], idxs[i,j+1,k]],           # left
        [idxs[i,j+1,k], idxs[i,j+1,k+1], idxs[i+1,j+1,k+1], idxs[i+1,j+1,k]]    # back
    ]
    
    add_volume!(_face_ids, topo)
end

# =============================================================================
# Test topology queries
# =============================================================================

println("Area node IDs for volume 1: ", get_area_node_ids(topo)[get_volume_area_ids(topo, 1)])
println("Volume node IDs for volume 2: ", get_volume_node_ids(topo, 2))

# =============================================================================
# Refine first volume and export
# =============================================================================
let topo = deepcopy(topo)
    @time for vol in RootIterator{3,4}(topo)
        if rand(0:1) |> Bool
            _refine!(vol,topo)
        end
    end
end

@time for vol in RootIterator{3,4}(topo)
    if rand(0:1) |> Bool
        _refine!(vol,topo)
    end
end

# @code_warntype _refine!(get_volumes(topo)[1],topo)



# vol = get_volumes(topo)[1]
# _refine!(vol, topo)
# vol = get_volumes(topo)[end]
# _refine!(vol, topo)

# _coarsen!(get_volumes(topo)[3],topo)

geometry_to_vtk(topo, "3d_test1")

# Attempt tetrahedralization for a random volume; export tets if successful
begin
    try
        vol_id = 1
        tets_local = tetrahedralize_volume_local_ids(topo, vol_id)
        tet_topo = build_tet_topology_from_volume(topo, vol_id; tets_local)
        geometry_to_vtk(tet_topo, "arrow_tets")
    catch err
        @warn "Ear-peeling tetrahedralization failed for 3d_test1 volume 1" err
    end
end

nothing