using OrderedCollections
using Ju3VEM
using Ju3VEM:FlattenVecs
using Chairmarks
using Random

@kwdef mutable struct ElementNodeMapping{DT<:AbstractDict}
    const dict::DT            = OrderedDict{Int,Int}()
    node_count          ::Int = 0
    edge_node_count     ::Int = 0
    face_moment_count   ::Int = 0
    volume_moment_count ::Int = 0
    isertion_finished   ::BitVector = falses(4)
end

function get_local_id(vnm::ElementNodeMapping) 
    return  vnm.node_count          + 
            vnm.edge_node_count     + 
            vnm.face_moment_count   + 
            vnm.volume_moment_count
end
function get_n_nodes(vnm::ElementNodeMapping) 
    return get_local_id(vnm)
end


# Generic function
function add_local_id!(vnm::ElementNodeMapping, node_id::Int, ::Val{N}) where N
    if haskey(vnm.dict, node_id) 
        return vnm.dict[node_id]
    end

    N > 1 &&!vnm.isertion_finished[N-1] && error("Insertion not finished for previous type $(N-1)")
    vnm.isertion_finished[N] == true    && error("Insertion already finished for type $N")

    _increment_counter!(vnm, Val(N))
    local_id = get_local_id(vnm)
    vnm.dict[node_id] = local_id

    return local_id
end

function get_local_id!(vnm::ElementNodeMapping, node_id::Int)
    vnm.dict[node_id]
end

# Dispatch for incrementing the right counter
_increment_counter!(vnm::ElementNodeMapping, ::Val{1}) = vnm.node_count += 1
_increment_counter!(vnm::ElementNodeMapping, ::Val{2}) = vnm.edge_node_count += 1
_increment_counter!(vnm::ElementNodeMapping, ::Val{3}) = vnm.face_moment_count += 1
_increment_counter!(vnm::ElementNodeMapping, ::Val{4}) = vnm.volume_moment_count += 1

# Convenience functions with original names
add_node_id!(vnm::ElementNodeMapping, node_id::Int) = 
    add_local_id!(vnm, node_id, Val(1))
add_edge_node_id!(vnm::ElementNodeMapping, node_id::Int) = 
    add_local_id!(vnm, node_id, Val(2))
add_face_moment_id!(vnm::ElementNodeMapping, node_id::Int) = 
    add_local_id!(vnm, node_id, Val(3))
add_volume_moment_id!(vnm::ElementNodeMapping, node_id::Int) = 
    add_local_id!(vnm, node_id, Val(4))


@inline finish_insertion!(vnm::ElementNodeMapping, 
            ::Val{N}) where N = vnm.isertion_finished[N] = true

@inline finish_node_insertion!(vnm::ElementNodeMapping) = finish_insertion!(vnm, Val(1))
@inline finish_edge_node_insertion!(vnm::ElementNodeMapping) = finish_insertion!(vnm, Val(2))
@inline finish_face_moment_insertion!(vnm::ElementNodeMapping) = finish_insertion!(vnm, Val(3))
@inline finish_volume_moment_insertion!(vnm::ElementNodeMapping) = finish_insertion!(vnm, Val(4))


function get_node_ids_view(enm::ElementNodeMapping)
    return view(enm.dict.keys,1:enm.node_count) 
end

function get_local_node_ids_view(enm::ElementNodeMapping)
    return view(enm.dict.vals,1:enm.node_count) 
end

function get_edge_node_ids_view(enm::ElementNodeMapping)
    start_idx = enm.node_count + 1
    end_idx = start_idx + enm.edge_node_count - 1
    return view(enm.dict.keys,start_idx:end_idx)
end

function get_local_edge_node_ids_view(enm::ElementNodeMapping)
    start_idx = enm.node_count + 1
    end_idx = start_idx + enm.edge_node_count - 1
    return view(enm.dict.vals,start_idx:end_idx)
end

function get_face_moment_ids_view(enm::ElementNodeMapping)
    start_idx = enm.node_count + enm.edge_node_count + 1
    end_idx = start_idx + enm.face_moment_count - 1
    return view(enm.dict.keys,start_idx:end_idx)
end

function get_local_face_moment_ids_view(enm::ElementNodeMapping)
    start_idx = enm.node_count + enm.edge_node_count + 1
    end_idx = start_idx + enm.face_moment_count - 1
    return view(enm.dict.vals,start_idx:end_idx)
end

function get_volume_moment_ids_view(enm::ElementNodeMapping)
    start_idx = enm.node_count + enm.edge_node_count + enm.face_moment_count + 1
    end_idx = start_idx + enm.volume_moment_count - 1
    return view(enm.dict.keys,start_idx:end_idx)
end

function get_local_volume_moment_ids_view(enm::ElementNodeMapping)
    start_idx = enm.node_count + enm.edge_node_count + enm.face_moment_count + 1
    end_idx = start_idx + enm.volume_moment_count - 1
    return view(enm.dict.vals,start_idx:end_idx)
end


# ids = collect(1:10)
# enm = ElementNodeMapping() 
# for id in ids
#     add_node_id!(enm,id)
# end
# finish_node_insertion!(enm)
# ids = collect(11:20)
# for id in ids
#     add_edge_node_id!(enm,id)
# end

# ids = collect(21:30)
# for id in ids
#     add_node_id!(enm,id)
# end









