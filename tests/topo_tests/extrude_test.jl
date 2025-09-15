using Ju3VEM
using StaticArrays
using Ju3VEM.VEMUtils: extrude_to_3d, create_voronoi_mesh
using Test

n = 10
nx,ny = n,n
nz = n
dx = 1/nx; dy = 1/ny

x_coords = range(0, nx*dx, length=nx+1)
y_coords = range(0, ny*dy, length=ny+1)

coords = [SA[x,y] for x in x_coords, y in y_coords]

topo = Topology{2}()
add_node!.(coords, Ref(topo))

idxs = LinearIndices((nx+1, ny+1))

for I in CartesianIndices((nx, ny))
    i, j = Tuple(I)
    node_ids = [idxs[i,j],idxs[i+1,j],idxs[i+1,j+1],idxs[i,j+1]]
    add_area!(node_ids,topo)
end

mesh = Mesh(topo, StandardEl{1}())

mesh3d = extrude_to_3d(nz,mesh)

@show length(mesh3d.topo.nodes)
@show length(unique(mesh3d.topo.nodes))

@test allunique(mesh3d.topo.nodes)

geometry_to_vtk(mesh3d.topo, "vtk/extrude_test")




mesh2d = create_voronoi_mesh((0.0,0.0),(1.0,1.0),nx,ny,StandardEl{1})

mesh3d = extrude_to_3d(nz,mesh2d)

@test allunique(mesh3d.topo.nodes)

geometry_to_vtk(mesh3d.topo, "vtk/extrude_test_voronoi")