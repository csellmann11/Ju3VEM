using StaticArrays
using ..VEMGeo: FaceData
using OrderedCollections: OrderedDict

struct NodeID2LocalID 
    map                  ::OrderedDict{Int32,Int16}
    counter              ::MVector{4,Int16} 
end

function NodeID2LocalID(;sizehint::Integer=100) 
    map = Dict{Int32,Int16}()
    sizehint!(map,sizehint)
    NodeID2LocalID(map,zero(MVector{4,Int16}))
end

@inline get_n_nodes(ntl::NodeID2LocalID) = length(ntl.map)

function add_field!(ntl::NodeID2LocalID,
    key::Integer,
    ::Val{N}) where N

    id = get!(ntl.map,key) do 
        ntl.counter[N]   += 1
        get_n_nodes(ntl) + 1
    end
    id
end

@inline add_vertex_id!(ntl::NodeID2LocalID,key::Integer) = add_field!(ntl,key,Val(1))
@inline add_edge_vertex_id!(ntl::NodeID2LocalID,key::Integer) = add_field!(ntl,key,Val(2))
@inline add_face_moment_id!(ntl::NodeID2LocalID,key::Integer) = add_field!(ntl,key,Val(3))
@inline add_volume_moment_id!(ntl::NodeID2LocalID,key::Integer) = add_field!(ntl,key,Val(4))


@inline get_local_id(ntl::NodeID2LocalID,key::Integer) = ntl.map[key]


# TODO: remove dependency on type FaceData
function create_node_mapping!(
    ntl::NodeID2LocalID,
    volume_id::Integer,mesh::Mesh{D,ET},
    facedata_col::Dict{<:Integer,<:FaceData{D,L}}) where {D,L,K,ET<:ElType{K}}


    empty!(ntl.map)
    ntl.counter .= 0

    topo = mesh.topo 
    
    iterate_volume_areas(facedata_col,topo,volume_id) do face, face_data, _
        node_ids = face_data.face_node_ids
        for id in node_ids.v.args[1]
            add_vertex_id!(ntl,id)
        end
        K == 1 && return 
        for id in node_ids.v.args[2]
            add_edge_vertex_id!(ntl,id)
        end

        for id in node_ids.v.args[3]
            add_face_moment_id!(ntl,id)
        end
    end

    for id in mesh.int_coords_connect[4][volume_id]
        add_volume_moment_id!(ntl,id)
    end
    return ntl
end

function create_node_mapping(volume_id,mesh::Mesh,
    facedata_col::Dict) 

    ntl = NodeID2LocalID(;sizehint = 30) 

    create_node_mapping!(ntl,volume_id,mesh,facedata_col)
    return ntl
end




