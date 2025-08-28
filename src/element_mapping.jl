# Notes for refactoring: 
# Dependencies: 
# - FaceData
# - Mesh



using StaticArrays


struct NodeID2LocalID 
    map                  ::Dict{Int,Int16}
    counter              ::MVector{4,Int16} 
end

function NodeID2LocalID(;sizehint::Int=100) 
    map = Dict{Int,Int16}()
    sizehint!(map,sizehint)
    NodeID2LocalID(map,zero(MVector{4,Int16}))
end

@inline get_n_nodes(ntl::NodeID2LocalID) = length(ntl.map)

function add_field!(ntl::NodeID2LocalID,
    key::Int,
    ::Val{N}) where N

    id = get!(ntl.map,key) do 
        ntl.counter[N]   += 1
        get_n_nodes(ntl) + 1
    end
    id
end

@inline add_vertex_id!(ntl::NodeID2LocalID,key::Int) = add_field!(ntl,key,Val(1))
@inline add_edge_vertex_id!(ntl::NodeID2LocalID,key::Int) = add_field!(ntl,key,Val(2))
@inline add_face_moment_id!(ntl::NodeID2LocalID,key::Int) = add_field!(ntl,key,Val(3))
@inline add_volume_moment_id!(ntl::NodeID2LocalID,key::Int) = add_field!(ntl,key,Val(4))


@inline get_local_id(ntl::NodeID2LocalID,key::Int) = ntl.map[key]



function create_node_mapping(volume_id::Int,mesh::Mesh{D,ET},
    facedata_col::Dict{Int,FaceData{D,L}}) where {D,L,K,ET<:ElType{K}}

    topo = mesh.topo 
    ntl = NodeID2LocalID(;sizehint = 30) 
    iterate_volume_areas(facedata_col,topo,volume_id) do _, face_data, _
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
  
        D == 2 && return 
        for id in node_ids.v.args[4]
            add_volume_moment_id!(ntl,id)
        end
    end
    return ntl
end







