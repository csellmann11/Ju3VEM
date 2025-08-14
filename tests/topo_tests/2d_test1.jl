using StaticArrays, WriteVTK
using OrderedCollections, Bumper
using LinearAlgebra, Statistics
using SmallCollections, Chairmarks



# Include the main module (professional approach)
include("../../src/small_utils.jl")
include("../../src/topo.jl")
include("../../src/vtkexports.jl")
include("../../src/element_refinement.jl")
include("../../src/lagrange_utils.jl")
include("../../src/element_coarsening.jl")
include("../../src/triangulation.jl")
# =============================================================================
# 3D Topology Test - Experimental Code
# =============================================================================

# Grid parameters
nx =30*30; ny = 30
dx = 1.0; dy = 1.0

# Generate coordinate ranges
x_coords = range(0, nx*dx, length=nx+1)
y_coords = range(0, ny*dy, length=ny+1)

# Create coordinate points
coords = [SA[x,y] for x in x_coords, y in y_coords]

# Initialize topology
topo = Topology{2}()

# Add all nodes to topology
add_node!.(coords, Ref(topo))

# Create linear indices for easy access
idxs = LinearIndices((nx+1, ny+1))

# Test coordinate access
println("Coordinate at (2,2): ", coords[idxs[2,2]])

# =============================================================================
# Create hexahedral elements
# =============================================================================

for I in CartesianIndices((nx, ny))
    i, j = Tuple(I)
    
    # Define face IDs for each hexahedron
    node_ids = [idxs[i,j],idxs[i+1,j],idxs[i+1,j+1],idxs[i,j+1]]
    add_area!(node_ids,topo)
end

# for element in RootIterator{2,3}(topo)
#     println(element)
# end

let topo = deepcopy(topo)
    @time for vol in RootIterator{2,3}(topo)
        if rand(0:1) |> Bool
            _refine!(vol,topo)
        end
    end
end



# =============================================================================
# Refine the topology
# =============================================================================
# element = get_areas(topo)[1]
# @time _refine!(element,topo)
# element = get_areas(topo)[end]
# _refine!(element,topo)
# element = get_areas(topo)[end]
# _refine!(element,topo)

# element = get_areas(topo)[2]
# _refine!(element,topo)

# _coarsen!(get_areas(topo)[3],topo)

# Export to VTK for visualization
# =============================================================================

geometry_to_vtk(topo, "2d_test1")

