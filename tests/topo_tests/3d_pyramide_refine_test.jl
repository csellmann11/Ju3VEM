using StaticArrays, WriteVTK
using OrderedCollections, Bumper
using LinearAlgebra, Statistics
using SmallCollections, Chairmarks
using Ju3VEM

# Build a square pyramid: base in z=0 plane, apex above center
topo = Topology{3}()

# Base square corners
base = [SA[0.0, 0.0, 0.0],
        SA[1.0, 0.0, 0.0],
        SA[1.0, 1.0, 0.0],
        SA[0.0, 1.0, 0.0]]
apex = SA[0.5, 0.5, 1.0]

node_ids = Int[]
append!(node_ids, add_node!.(base, Ref(topo)))
push!(node_ids, add_node!(apex, topo))

# Faces: base quad (ccw) and four triangular sides
base_ids = node_ids[1:4]
apex_id = node_ids[5]

faces = [
    base_ids,
    [base_ids[1], base_ids[2], apex_id],
    [base_ids[2], base_ids[3], apex_id],
    [base_ids[3], base_ids[4], apex_id],
    [base_ids[4], base_ids[1], apex_id]
]

add_volume!(faces, topo)
 


# Refine the pyramid volume 


vol = get_volumes(topo)[1]
_refine!(vol, topo)

vol = get_volumes(topo)[end]
_coarsen!(vol,topo)
_coarsen!(vol,topo)
# _coarsen!(get_volumes(topo)[5],topo)

# vol = get_volumes(topo)[1]
# _refine!(vol, topo)



geometry_to_vtk(topo, "vtk/3d_pyramide_refine_test")

# Also test tetrahedralization of the last volume if possible
try
    vol_id = length(get_volumes(topo))
    tets_local, local_to_global = tetrahedralize_volume(topo, vol_id)
    tet_topo = build_tet_topology_from_volume(topo, vol_id; tets_local)
    geometry_to_vtk(tet_topo, "pyramid_tets")
catch err
    @warn "Ear-peeling tetrahedralization failed for pyramid test (acceptable for non-star-shaped cases)" err
end


