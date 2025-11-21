using StaticArrays, WriteVTK
using OrderedCollections, Bumper
using LinearAlgebra, Statistics
using SmallCollections, Chairmarks
using LoopVectorization
using Random
using Ju3VEM
using Ju3VEM.VEMGeo: _refine!,_coarsen!,refine!
# =============================================================================
# 3D Topology Test - Experimental Code
# =============================================================================

# Grid parameters
# nx = 30; ny = 30; nz = 30
nx,ny,nz = 30,30,30
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
ix(i, j, k) = Int32(idxs[i, j, k])
for I in CartesianIndices((nx, ny, nz))
    i, j, k = Tuple(I)
    
    # Define face IDs for each hexahedron
    _face_ids = [
        [ix(i,j,k), ix(i+1,j,k), ix(i+1,j+1,k), ix(i,j+1,k)],           # bottom
        [ix(i,j,k+1), ix(i+1,j,k+1), ix(i+1,j+1,k+1), ix(i,j+1,k+1)],   # top
        [ix(i,j,k), ix(i+1,j,k), ix(i+1,j,k+1), ix(i,j,k+1)],           # front 
        [ix(i+1,j,k), ix(i+1,j+1,k), ix(i+1,j+1,k+1), ix(i+1,j,k+1)],   # right 
        [ix(i,j,k), ix(i,j,k+1), ix(i,j+1,k+1), ix(i,j+1,k)],           # left
        [ix(i,j+1,k), ix(i,j+1,k+1), ix(i+1,j+1,k+1), ix(i+1,j+1,k)]    # back
    ]
    
    add_volume!(_face_ids, topo)
end




# =============================================================================
# Test topology queries
# =============================================================================

println("Area node IDs for volume 1: ", get_area_node_ids(topo)[get_volume_area_ids(topo, 1)])
println("Volume node IDs for volume 2: ", get_volume_node_ids(topo, 2))


# Ju3VEM.VEMGeo.get_volume_node_ids(topo, 1)
# =============================================================================
# Refine first volume and export
# =============================================================================
using JET
let topo = deepcopy(topo)
    GC.gc(true)
    rng = Random.MersenneTwister(123)
    @time for vol in RootIterator{3,4}(topo)
        if rand(rng,0:1) |> Bool
            _refine!(vol,topo)
        end
    end


    # @time for vol in RootIterator{3,4}(topo)
    #     if rand(rng,0:1) |> Bool
    #         _refine!(vol,topo)
    #     end
    # end
    # # test of refinement warpper 
    elements_to_refine = BitVector(rand(rng,0:1) |> Bool for _ in get_volumes(topo))
    els_to_ref = sum(elements_to_refine)
    println("Elements to refine: $els_to_ref")
    # @time refine!(topo, elements_to_refine)
    # t = @report_opt refine!(topo, elements_to_refine)
    # show(t)

    @show length(RootIterator{4}(topo))
    write_vtk(topo, "vtk/3d_test1")
end




# rng = Random.MersenneTwister(123)
# elements_to_refine = BitVector(rand(rng,0:1) |> Bool for _ in get_volumes(topo))
# @code_warntype create_refined_topology_holistic(topo, elements_to_refine)


# @code_warntype _refine!(get_volumes(topo)[1],topo)



# vol = get_volumes(topo)[1]
# _refine!(vol, topo)
# vol = get_volumes(topo)[end]
# _refine!(vol, topo)

# _coarsen!(get_volumes(topo)[3],topo)

# write_vtk(topo, "vtk/3d_test1")

