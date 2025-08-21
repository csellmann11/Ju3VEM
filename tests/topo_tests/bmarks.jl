using Chairmarks
using StaticArrays, WriteVTK
using OrderedCollections, Bumper
using LinearAlgebra, Statistics

using Ju3VEM




# Build a square pyramid: base in z=0 plane, apex above center

function create_pyramid()
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
    vol = get_volumes(topo)[1]
    return vol,topo
end

# Refine the pyramid volume




# _refine!(vol, topo)
# find_single_intersec(get_volume_node_ids(topo),vol.childs)

# @b find_single_intersec(get_volume_node_ids($topo),$vol.childs)

# a = ones(4,4)

# @b create_pyramid() _refine!(_...)

# vol,topo = create_pyramid();
# @time _refine!(vol,topo)


# v = [rand(1:20,4) for _ in 1:8]

# @b get_unique_values($v,1:length($v))


