

struct Node{D,T} <: StaticVector{D,T}
    id::Int32
    coords::SVector{D,T}
end

function Node(id::Integer, coords::StaticVector{D,T}) where {D,T}
    Node{D,T}(Int32(id), coords)
end

Node(coords::AbstractVector) = Node(-1, coords)
Node(coords::StaticVector) = Node(-1, coords)



get_coords(n::Node) = n.coords
get_id(n::Node) = n.id
get_id(i::Integer) = i
Base.getindex(n::Node{D,T}, i::Int32) where {D,T} = get_coords(n)[i]
Base.getindex(n::Node{D,T}, i::Int64) where {D,T} = get_coords(n)[i]

Base.zero(::Type{Node{D,T}}) where {D,T} = Node(zero(SVector{D,T})) # to create placeholder nodes
Base.zero(n::Node) = zero(typeof(n))
Base.rand(::Type{Node{D,T}}) where {D,T} = Node(-1, @SVector rand(T,D))


#########################
# NManifold
#########################
@kwdef struct NManifold{N,D}
    id::Int32
    refinement_level::Int32

    # tree structure
    parent_id::Int32
    childs::Vector{Int32} = Int32[]
end

const Edge{D} = NManifold{2,D}
const Area{D} = NManifold{3,D}
const Volume{D} = NManifold{4,D}

function NManifold{N,D}(id::Integer, refinement_level::Integer, parent_id::Integer) where {N,D}
    NManifold{N,D}(id=Int32(id), refinement_level=Int32(refinement_level), parent_id=Int32(parent_id))
end

# function activate!(mf::MF) where MF <: Union{<:NManifold,<:Node}
#     mf.id = abs(mf.id)
#     return mf
# end



# is_active_root(mf::NManifold) = is_root(mf) && is_active(mf)

#########################  
# Topology
#########################
using FixedSizeArrays
Base.convert(::FixedSizeVector{T,V},v::SVector{N,T}) where {T,N,V} = FixedSizeVector{T,V}(v) 
const ID_VEC_TYPE{T} = FixedSizeVector{T,Memory{T}}
@kwdef struct Topology{D}

    nodes::Vector{Node{D,Float64}}              = Node{D,Float64}[]
    edges::OrderedDict{ID_VEC_TYPE{Int32},Edge{D}}     = OrderedDict{ID_VEC_TYPE{Int32},Edge{D}}()
    areas::OrderedDict{ID_VEC_TYPE{Int32},Area{D}}     = OrderedDict{ID_VEC_TYPE{Int32},Area{D}}()
    volumes::OrderedDict{ID_VEC_TYPE{Int32},Volume{D}} = OrderedDict{ID_VEC_TYPE{Int32},Volume{D}}()

    is_active_node::Vector{Bool} = Bool[]

    is_root_edge::Vector{Bool} = Bool[]
    is_active_edge::Vector{Bool} = Bool[]
    
    is_root_area::Vector{Bool} = Bool[]
    is_active_area::Vector{Bool} = Bool[]

    is_root_volume::Vector{Bool} = Bool[]
    is_active_volume::Vector{Bool} = Bool[]

    connectivity::Matrix{Vector{ID_VEC_TYPE{Int32}}} = begin
        mat = Matrix{Vector{ID_VEC_TYPE{Int32}}}(undef, 4, 4)
        mat[1, 2] = edges.keys
        mat[1, 3] = areas.keys
        mat[3, 4] = volumes.keys #! volumes are index by face_ids
        mat[1, 4] = Vector{ID_VEC_TYPE{Int32}}()
        for i in 1:D+1, j in 2:i-1
            if i != j && (i,j) != (4,3)
                mat[j, i] = Vector{ID_VEC_TYPE{Int32}}()
            end
        end 
        mat
    end
end

@inline get_nodes(topo::Topology{D}) where D = topo.nodes::Vector{Node{D,Float64}}
@inline get_edges(topo::Topology{D}) where D = topo.edges.vals::Vector{Edge{D}}
@inline get_areas(topo::Topology{D}) where D = topo.areas.vals::Vector{Area{D}}
@inline get_volumes(topo::Topology{D}) where D = topo.volumes.vals::Vector{Volume{D}}


@inline get_geo_vec(topo, ::Val{1}) = get_nodes(topo)
@inline get_geo_vec(topo, ::Val{2}) = get_edges(topo)
@inline get_geo_vec(topo, ::Val{3}) = get_areas(topo)
@inline get_geo_vec(topo, ::Val{4}) = get_volumes(topo)
function get_geo_vec(::Any, ::Val{N}) where {N}
    error("N must be 1,2,3 or 4")
end

#Info: access volume information
@inline get_volume_node_ids(topo)::Vector{ID_VEC_TYPE{Int32}}      = topo.connectivity[1, 4]
@inline get_volume_node_ids(topo, id::Integer)::ID_VEC_TYPE{Int32}     = get_volume_node_ids(topo)[id]
@inline get_volume_edge_ids(topo)::Vector{ID_VEC_TYPE{Int32}}      = topo.connectivity[2, 4]
@inline get_volume_edge_ids(topo, id::Integer)::ID_VEC_TYPE{Int32}     = get_volume_edge_ids(topo)[id]
@inline get_volume_area_ids(topo)::Vector{ID_VEC_TYPE{Int32}}      = topo.connectivity[3, 4]
@inline get_volume_area_ids(topo, id::Integer)::ID_VEC_TYPE{Int32}     = get_volume_area_ids(topo)[id]
#Info: access face information
@inline get_area_node_ids(topo)::Vector{ID_VEC_TYPE{Int32}}        = topo.connectivity[1, 3]
@inline get_area_node_ids(topo, id::Integer)::ID_VEC_TYPE{Int32}       = get_area_node_ids(topo)[id]
@inline get_area_edge_ids(topo)::Vector{ID_VEC_TYPE{Int32}}        = topo.connectivity[2, 3]
@inline get_area_edge_ids(topo, id::Integer)::ID_VEC_TYPE{Int32}       = get_area_edge_ids(topo)[id]
#Info: access edge information
@inline get_edge_node_ids(topo)::Vector{ID_VEC_TYPE{Int32}}        = topo.connectivity[1, 2]
@inline get_edge_node_ids(topo, id::Integer)::ID_VEC_TYPE{Int32}       = get_edge_node_ids(topo)[id]


@inline add_volume_node_id!(topo, dest::Integer, val::Integer) =
    push!(get_volume_node_ids(topo, dest), val)
@inline add_volume_edge_id!(topo, dest::Integer, val::Integer) =
    push!(get_volume_edge_ids(topo, dest), val)
@inline add_volume_area_id!(topo, dest::Integer, val::Integer) =
    push!(get_volume_area_ids(topo, dest), val)

@inline add_area_node_id!(topo, dest::Integer, val::Integer) =
    push!(get_area_node_ids(topo, dest), val)
@inline add_area_node_id!(topo,loc::Integer,dest::Integer,val::Integer) = 
    get_area_node_ids(topo,dest)[loc] = val


@inline add_area_edge_id!(topo, dest::Integer, val::Integer) =
    push!(get_area_edge_ids(topo, dest), val)
@inline add_area_edge_id!(topo,loc::Integer,dest::Integer,val::Integer) = 
    get_area_edge_ids(topo,dest)[loc] = val


@inline add_edge_node_id!(topo, dest::Integer, val::Integer) =
    push!(get_edge_node_ids(topo, dest), val)
@inline add_edge_node_id!(topo,loc::Integer,dest::Integer,val::Integer) = 
    get_edge_node_ids(topo,dest)[loc] = val


@inline __default_node_id(topo) = length(get_nodes(topo)) + 1



function change_activity!(v::AbstractVector{Bool},id::Integer,activity::Bool)
    v[id] = activity
    return nothing
end

@inline change_activity!(n::Node,topo::Topology,activity::Bool) = change_activity!(topo.is_active_node,n.id,activity)
@inline change_activity!(m::NManifold{2},topo::Topology,activity::Bool) = change_activity!(topo.is_active_edge,m.id,activity)
@inline change_activity!(m::NManifold{3},topo::Topology,activity::Bool) = change_activity!(topo.is_active_area,m.id,activity)
@inline change_activity!(m::NManifold{4},topo::Topology,activity::Bool) = change_activity!(topo.is_active_volume,m.id,activity)

@inline activate!(m,topo::Topology) = change_activity!(m,topo,true)
@inline deactivate!(m,topo::Topology) = change_activity!(m,topo,false)



function change_root_status!(v::AbstractVector{Bool},id::Integer,status::Bool)
    v[id] = status
    return nothing
end

@inline change_root_status!(m::NManifold{2},topo::Topology,status::Bool) = change_root_status!(topo.is_root_edge,m.id,status)
@inline change_root_status!(m::NManifold{3},topo::Topology,status::Bool) = change_root_status!(topo.is_root_area,m.id,status)
@inline change_root_status!(m::NManifold{4},topo::Topology,status::Bool) = change_root_status!(topo.is_root_volume,m.id,status)

@inline make_root!(m,topo::Topology) = change_root_status!(m,topo,true)
@inline remove_root!(m,topo::Topology) = change_root_status!(m,topo,false)

function init_activity!(topo::Topology,id::Integer)
    if id > length(topo.is_active_node)
        resize_with_trailing_zeros!(topo.is_active_node, 2*id)
    end
    topo.is_active_node[id] = true
    return nothing
end

function init_activity_and_root!(topo::Topology,id::Integer,::Val{N}) where N
    root_vec,activity_vec = if N == 2
        topo.is_root_edge,topo.is_active_edge
    elseif N == 3
        topo.is_root_area,topo.is_active_area
    elseif N == 4
        topo.is_root_volume,topo.is_active_volume
    end
    if id > length(root_vec)
        resize_with_trailing_zeros!(root_vec, 2*id)
    end
    root_vec[id] = true

    if id > length(activity_vec)
        resize_with_trailing_zeros!(activity_vec, 2*id)
    end
    
    activity_vec[id] = true
    return nothing
end
    






get_id(mf::NManifold) = mf.id
@inline is_active(n::Node,topo::Topology) = topo.is_active_node[n.id]
@inline is_active(m::NManifold{2},topo::Topology) = topo.is_active_edge[m.id]
@inline is_active(m::NManifold{3},topo::Topology) = topo.is_active_area[m.id]
@inline is_active(m::NManifold{4},topo::Topology) = topo.is_active_volume[m.id]







function add_node!(coords::AbstractVector,
    topo::Topology{D},
    id=__default_node_id(topo)) where {D}

    node = Node(id, coords)
    nodes = get_nodes(topo)
    push!(nodes, node)
    init_activity!(topo,id)
    return id
end
function add_node!(node::Node,
    topo::Topology)
    add_node!(node.coords, topo, node.id)

    return node.id
end


function rotate_to_min(v::AbstractVector{<:Integer})
    # Rotate the vector so that the smallest element is first
    dest = ID_VEC_TYPE{Int32}(undef, length(v))
    circshift!(dest,v,-(argmin(v) - 1))
    return dest
end

@inline edge_hash(n1::Integer, n2::Integer) = ID_VEC_TYPE(SVector(minmax(Int32(n1), Int32(n2))))
@inline edge_hash(v::Union{SVector{2,I},Tuple{I,I}}) where I<:Integer = edge_hash(v...)
@inline area_hash(v) = rotate_to_min(v)
# @inline volume_hash(v) = sort!(ID_VEC_TYPE(v))
@inline volume_hash(v) = rotate_to_min(v)

function add_edge!(node_ids::AbstractVector{<:Integer},
    topo::Topology{D},
    parent_id::Integer=Int32(0),
    parent_ref_level::Integer=Int32(0)) where D

    @assert length(node_ids) == 2

    temp_id = length(topo.edges) + 1
    # edge = get!(topo.edges, edge_hash(node_ids),
    #     Edge{D}(temp_id, parent_ref_level + 1, parent_id))

    generated_hash = edge_hash(node_ids)
    idx = OrderedCollections.ht_keyindex2(topo.edges, generated_hash)
    id = if idx <= 0
        topo.edges[generated_hash] = Edge{D}(temp_id, parent_ref_level + 1, parent_id)
        init_activity_and_root!(topo,temp_id,Val(2))
        temp_id
    else 
        edge = topo.edges.vals[idx]
        activate!(edge,topo)
        edge.id
    end
    
    return id
end

function full_activate!(area::Area{D}, topo::Topology{D})::Nothing where D
    area.id > 0 && return nothing
    activate!(area,topo)
    for edge_id in get_area_edge_ids(topo,area.id)
        edge = get_edges(topo)[edge_id]
        activate!(edge,topo)
    end
    nothing
end

function add_area!(_node_ids::AbstractVector{<:Integer},
    topo::Topology{D},
    parent_id::Integer=Int32(0),
    parent_ref_level::Integer=Int32(0)) where D

    area_dict = topo.areas
    temp_id = length(area_dict) + 1


    # Try to find area by hash (normal or reversed)
    generated_hash = area_hash(_node_ids)
    idx = OrderedCollections.ht_keyindex2(area_dict, generated_hash)
    
    if idx <= 0
        # Not found with normal hash, try reversed
        rev_hash = area_hash(reverse(_node_ids))
        idx = OrderedCollections.ht_keyindex2(area_dict, rev_hash)
    end
    
    # Create new area if not found
    area_found = (idx > 0)
    if !area_found
        area_dict[generated_hash] = Area{D}(temp_id, parent_ref_level + 1, parent_id)
        init_activity_and_root!(topo,temp_id,Val(3))
        idx = temp_id
    end
    
    area = area_dict.vals[idx]

 
    !is_active(area,topo) && full_activate!(area, topo)
    id = area.id
    # id != temp_id && return area.id
    area_found && return area.id

    node_ids = get_area_node_ids(topo, id)


    area_edge_ids = ID_VEC_TYPE{Int32}(undef,length(node_ids))
    push!(get_area_edge_ids(topo), area_edge_ids)
    for (i, node_id) in enumerate(node_ids)
        ip1 = get_next_idx(node_ids, i)
        edge_id = add_edge!(SA[node_id, node_ids[ip1]],
            topo, Int32(0) , parent_ref_level)

        area_edge_ids[i] = edge_id
    end
    return id
end

function full_activate!(volume::Volume{D}, topo::Topology{D})::Nothing where D
    volume.id > 0 && return nothing
    activate!(volume,topo)

    for area_id in get_volume_area_ids(topo,volume.id)
        area = get_areas(topo)[area_id]
        activate!(area,topo)
    end
    nothing
end


function add_volume!(face_ids::AbstractVector{<:Integer},
    topo::Topology{D},
    parent_id::Integer=Int32(0),
    parent_ref_level::Integer=Int32(0)) where D
    
    volume_dict = topo.volumes
    temp_id = length(volume_dict) + 1

    volume = get!(volume_dict, volume_hash(face_ids),
        Volume{D}(temp_id, parent_ref_level + 1, parent_id))
    
    if volume.id == temp_id
        init_activity_and_root!(topo,temp_id,Val(4))
    end

    !is_active(volume,topo) && full_activate!(volume, topo)
    id = volume.id 

    
    if id == temp_id
        # unique_ids = get_unique_values(get_area_node_ids(topo),face_ids)
        area_node_ids_full = get_area_node_ids(topo) 

        id_set = Set(Iterators.flatten(area_node_ids_full[face_ids]))
        ids = FixedSizeVector{Int32}(undef, length(id_set))
        copyto!(ids, id_set)

        push!(get_volume_node_ids(topo), ids)
    end
    volume.id

end


function add_volume!(node_ids_col::AbstractVector{<:AbstractVector{<:Integer}},
    topo::Topology{D},
    parent_id::Integer=Int32(0),
    parent_ref_level::Integer=Int32(0)) where D

    n_faces = length(node_ids_col)
    volume_dict = topo.volumes
    temp_id = length(volume_dict) + 1
 

    _face_ids = ID_VEC_TYPE{Int32}(undef, n_faces)

    for (i, node_ids) in enumerate(node_ids_col)
        face_id = add_area!(node_ids, topo, Int32(0), parent_ref_level)
        _face_ids[i] = face_id
    end

    volume = get!(volume_dict, volume_hash(_face_ids),
        Volume{D}(temp_id, parent_ref_level + 1, parent_id))

    if volume.id == temp_id
        init_activity_and_root!(topo,temp_id,Val(4))
    end

    activate!(volume, topo)
    id = volume.id


    if id == temp_id
        push!(get_volume_node_ids(topo), get_unique_values(node_ids_col))
    end
    volume.id
    
end


@inline is_root(m::NManifold{2},topo::Topology) = topo.is_root_edge[m.id]
@inline is_root(m::NManifold{3},topo::Topology) = topo.is_root_area[m.id]
@inline is_root(m::NManifold{4},topo::Topology) = topo.is_root_volume[m.id]

# @inline is_root(mf::NManifold) = mf.is_root

@inline is_active_root(mf::NManifold,topo::Topology) = is_root(mf,topo) && is_active(mf,topo)

#TODO: check performance implications of union return type
@inline function get_edge(n1::Int32, n2::Int32, topo::Topology{D})::Edge{D} where D
    edge = get(topo.edges, edge_hash(SVector(n1, n2)), nothing)
    @assert edge !== nothing "Edge not found"
    return edge
end

@inline function has_edge(n1::Int32, n2::Int32, topo::Topology{D}) where D
    edge_hash = edge_hash(SVector(n1, n2))
    return haskey(topo.edges, edge_hash)
end
@inline has_edge(n1::Node, n2::Node, topo::Topology) = has_edge(n1.id, n2.id, topo)

@inline get_edge(n1::Node, n2::Node, topo::Topology) = get_edge(n1.id, n2.id, topo)

@inline function get_face(v::AbstractVector{Int32}, topo::Topology{D}) where D
    area_hash = area_hash(v)
    return get_areas(topo)[area_hash]
end
@inline get_face(v::AbstractVector{<:Edge}, topo::Topology) = get_face(get_id.(v), topo)














# ======================================================================
# Iterator for the topology
# ======================================================================
"""
    RootIterator{D,Idx}(topo::Topology{D})

Iterator to iterate over all active roots in a topology.
Freezes the topology, so no new elements can be added or coarsed.

E.g. 
- `RootIterator{D,3}` Iterates over all Areas
- `RootIterator{D,2}` Iterates over all Edges
"""
struct RootIterator{D,Idx}
    topo::Topology{D}
    state_is_root::BitVector
end

function RootIterator{Idx}(topo::Topology{D}) where {D,Idx}
    return RootIterator{D,Idx}(topo)
end

function RootIterator{D,Idx}(topo::Topology{D}) where {D,Idx}
    state_is_root::BitVector = is_active_root.(get_geo_vec(topo, Val(Idx)),Ref(topo))
    return RootIterator{D,Idx}(topo, state_is_root)
end

function RootIterator{D,1}(topo::Topology{D}) where {D}
    state_is_root::BitVector = is_active.(get_geo_vec(topo, Val(1)),Ref(topo))
    # state_is_root = trues(length(get_nodes(topo)))
    return RootIterator{D,1}(topo, state_is_root)
end

Base.firstindex(r::RootIterator{D,Idx}) where {D,Idx} = 1

function Base.iterate(r::RootIterator{D,Idx}, state=1) where {D,Idx}

    while state <= length(r.state_is_root)
        geo = get_geo_vec(r.topo, Val(Idx))[state]

        # is_active check is needed if an element has been coarsed 

        # When an element is coarsed also the other childs are coarsed
        # which are yet to be iterated
        if r.state_is_root[state] && is_active(geo,r.topo)
            return geo, state + 1
        end
        state += 1
    end
    return nothing
end

function get_num_active_nodes(topo::Topology)
    return sum(topo.is_active_node)
end

function get_num_active_root_edges(topo::Topology)
    return sum(a*r for (a,r) in zip(topo.is_active_edge,topo.is_root_edge))
end

function get_num_active_root_areas(topo::Topology)
    return sum(a*r for (a,r) in zip(topo.is_active_area,topo.is_root_area))
end

function get_num_active_root_volumes(topo::Topology)
    return sum(a*r for (a,r) in zip(topo.is_active_volume,topo.is_root_volume))
end

function get_num_active_roots(geo_vec::AbstractVector{T},topo::Topology) where T<:NManifold
    return count(geo -> is_root(geo,topo) && is_active(geo,topo), geo_vec)
end

# Base.length(r::RootIterator{D,Idx}) where {D,Idx} = sum(r.state_is_root)
# Base.length(r::RootIterator{D,Idx}) where {D,Idx} = get_num_active_roots(get_geo_vec(r.topo, Val(Idx)),r.topo) 
function Base.length(r::RootIterator{D,Idx}) where {D,Idx}
    if Idx == 1
        return get_num_active_nodes(r.topo)
    elseif Idx == 2
        return get_num_active_root_edges(r.topo)
    elseif Idx == 3
        return get_num_active_root_areas(r.topo)
    elseif Idx == 4
        return get_num_active_root_volumes(r.topo)
    end
end


#############################################
# Tree utils
#############################################
@kwdef mutable struct CustStack{T,V<:AbstractVector{T}} <: AbstractVector{T}
    const stack::V
    last ::Int32 = Int32(0)
end


Base.length(aq::CustStack) = aq.last
Base.size(aq::CustStack) = (aq.last,)

function Base.push!(aq::CustStack{T},entry::T) where T
    aq.last += 1
    @boundscheck (aq.last <= length(aq.stack))
    aq.stack[aq.last] = entry 
    nothing
end


function pop_last!(aq::CustStack) 
    aq.last -= 1
    return aq.stack[aq.last + 1]
end

Base.@propagate_inbounds function Base.getindex(aq::CustStack,i::Integer) 
    @boundscheck (i <= aq.last)
    return aq.stack[i]
end

@inline Base.isempty(aq::CustStack) = aq.last == 0

@inline function stack_to_fixed_size_vector(aq::CustStack{T}) where T
    v = FixedSizeVector{T}(undef,length(aq))
    for i in eachindex(v,aq)
        v[i] = aq[i]
    end
    return v
end



"""
    get_iterative_area_vertex_ids(area::Area{D}, topo::Topology{D},
        abs_ref_level::Int32=typemax(Int32)) where D

Get the node ids of the mesh of an area.

# Arguments
- `area::Area{D}`: The area to get the mesh node ids of.
- `topo::Topology{D}`: The topology to get the mesh node ids from.
- `abs_ref_level::Int32`: Does not go below this refinement level in the recursive search
"""
function get_iterative_area_vertex_ids!(vertex_ids::AbstractVector{Int32},
    area_id::Int32,
    topo::Topology{D},
    abs_ref_level::Int32=typemax(Int32)) where D
    
    cond(edge) = edge.refinement_level == abs_ref_level || is_root(edge,topo)
    iterate_element_edges(topo, area_id, cond) do node_id, _, _
        push!(vertex_ids, node_id)
    end
    return nothing
end

function get_iterative_area_vertex_ids!(vertex_ids::AbstractVector{Int32},
    area::Area{D},
    topo::Topology{D},
    abs_ref_level::Int32=typemax(Int32)) where D
    
    get_iterative_area_vertex_ids!(vertex_ids, area.id, topo, abs_ref_level)
end


function get_iterative_area_vertex_ids(area_id::Integer, 
    topo::Topology{D},
    abs_ref_level::Int32=typemax(Int32)) where D 

 
    @no_escape begin
        _vertex_ids = CustStack(stack = @alloc(Int32,100)) 
        get_iterative_area_vertex_ids!(_vertex_ids, area_id, topo, abs_ref_level)
        stack_to_fixed_size_vector(_vertex_ids)
    end
end


@inline get_iterative_area_vertex_ids(area::Area, 
                    topo,abs_ref_level::Int32=typemax(Int32)) =       
                        get_iterative_area_vertex_ids(area.id, topo, abs_ref_level)

"""
iterate_element_edges(f::F1,topo::Topology{D},area_id::Int32, cond::F2 = is_root)

Iterate over the edges of an element and apply a function on the edges

# Arguments
- `f::Function`: The function to apply on the edges. f takes 
the `node_id`, `edge_id` and `parent_edge_id` as arguments
"""
# function iterate_element_edges(fun::F1, topo::Topology{D}, area_id::Int32, cond::F2=is_root) where {D,F1,F2}
#     edge_ids = topo.connectivity[2, 3][area_id]
#     edges = get_edges(topo)

#     nodes_of_edges = get_edge_node_ids(topo)
#     area_nodes = get_area_node_ids(topo, area_id)

#     for (count, parent_edge_id) in enumerate(edge_ids)
#         edge = edges[parent_edge_id]

#         n1 = area_nodes[count]
#         countp1 = get_next_idx(area_nodes, count)
#         n2 = area_nodes[countp1]

#         apply_f_on(cond, edge, edges, n2 < n1) do root_edge, rev_child_order
#             # idx = rev_child_order ? 1 : 2
#             idx = 1 + rev_child_order
#             edge_id = root_edge.id 
#             node_id = nodes_of_edges[edge_id][idx]
#             fun(node_id, edge_id, parent_edge_id)
#         end
#     end
# end


function iterate_element_edges(fun::F1, 
    topo::Topology{D}, 
    area_id::Int32, 
    cond::F2=x -> is_root(x,topo)) where {D,F1,F2}
    
    edge_ids = get_area_edge_ids(topo, area_id)
    edges = get_edges(topo)

    nodes_of_edges = get_edge_node_ids(topo)
    area_node_ids = get_area_node_ids(topo, area_id)

    @no_escape begin
        aq = CustStack(stack = @alloc(Int32,100))

        last_node_id = first(area_node_ids)
        for parent_edge_id in edge_ids
            push!(aq,parent_edge_id)
 
            while !isempty(aq)
                root_edge_id = pop_last!(aq)
                root_edge = edges[root_edge_id]
                n1guess = get_edge_node_ids(topo, root_edge_id)[1] 
                rev = n1guess != last_node_id

                if cond(root_edge)
                    node_id      = nodes_of_edges[root_edge_id][1+rev]
                    last_node_id = nodes_of_edges[root_edge_id][2-rev]

                    fun(node_id,root_edge_id,parent_edge_id)
                    continue
                end

                # TODO: CHeck if this must be ifelse 
         
                push!(aq,root_edge.childs[2-rev])
                push!(aq,root_edge.childs[rev+1])
            end
        end
    end
end










# """
#     apply_f_on(f::F1, cond::F2,
#         feature::T, features::AbstractVector{T},
#         rev_child_order::Bool) where {T,F1<:Function,F2<:Function}

# Apply a function on a feature and its children if the condition is true.

# # Arguments 
# - `f::Function`: The function to apply on the feature. f takes 
# the `feature`, `rev_child_order` as arguments
# - `cond::Function`: The condition to check if the feature should be applied.
# - `feature::T`: The feature to apply the function on.
# - `features::AbstractVector{T}`: The vector of features to apply the function on.
# - `rev_child_order::Bool`: If true, the child order is reversed.
# """
# function apply_f_on(f::F1, cond::F2,
#     feature::T, features::AbstractVector{T},
#     rev_child_order::Bool) where {T,F1<:Function,F2<:Function}
#     if cond(feature) && is_active(feature)
#         f(feature, rev_child_order)
#         return nothing
#     end

#     # If this feature has no children, terminate traversal safely
#     # if isempty(feature.childs)
#     #     return nothing
#     # end
#     # transform_fun = ifelse(rev_child_order,reverse,identity)
#     # for child_id in transform_fun(feature.childs)
#     #     apply_f_on(f, cond, features[child_id], features, false)
#     # end

#     child_id1 = feature.childs[1+rev_child_order]
#     apply_f_on(f, cond, features[child_id1], features, false)

#     child_id2 = feature.childs[2-rev_child_order]
#     apply_f_on(f, cond, features[child_id2], features, true)

#     nothing
# end

# apply_f_on_roots(f, feature, features, rev_child_order::Bool=false) =
#     apply_f_on(f, is_root, feature, features, rev_child_order)


# function apply_f_on_unordered(f::F1, 
#     cond::F2, 
#     feature::T, 
#     features::AbstractVector{T}) where {F1<:Function,F2<:Function,T<:NManifold}

#     if cond(feature) && is_active(feature)
#         f(feature)
#         return nothing
#     end
    
#     for child_id in feature.childs 
#         apply_f_on_unordered(f, cond, features[child_id], features)
#     end
#     nothing
# end

# apply_f_on_unordered_roots(f, feature, features) =
#     apply_f_on_unordered(f, is_root, feature, features)



"""
iterate_volume_areas(fun::F1, 
    facedata_col::Dict{Int32,FD}, 
    topo::Topology{D}, volume_id::Int32, 
    cond::F2=is_root) where {D,F1,F2,FD} 

Iterate over the areas of a volume and apply a function on the areas.

# Arguments
- `fun::Function`: The function to apply on the areas. fun takes 
the `root_area`, `facedata` and `topo` as arguments
- `facedata_col::Dict{Int32,FD}`: The dictionary of face data.
- `topo::Topology{D}`: The topology to get the areas from.
- `volume_id::Int32`: The id of the volume to iterate over.
- `cond::Function`: The condition to check if the area should be applied.
"""
function iterate_volume_areas(fun::F1, 
    facedata_col::Dict{Int32,FD},
    topo::Topology{D},
    volume_id::Int32,
    cond::F2=x -> is_root(x,topo)) where {D,F1,F2,FD}
    
    area_ids = get_volume_area_ids(topo, volume_id)
    areas = get_areas(topo)
    @no_escape begin
        aq = CustStack(stack = @alloc(Int32,1000))
        
        for area_id in area_ids
            push!(aq,area_id)
            while !isempty(aq)
                root_area_id = pop_last!(aq)
                root_area = areas[root_area_id]
                
                if cond(root_area)
                    fun(root_area,facedata_col[root_area_id],topo)
                    continue
                end

                for child_id in root_area.childs
                    push!(aq,child_id)
                end
            end
        end
    end
end


function iterate_volume_areas(fun::F1, 
    topo::Topology{D},
    volume_id::Int32,
    cond::F2=x -> is_root(x,topo)) where {D,F1,F2}
    
    area_ids = get_volume_area_ids(topo, volume_id)
    areas = get_areas(topo)
    @no_escape begin
        aq = CustStack(stack = @alloc(Int32,1000))
        
        for area_id in area_ids
            push!(aq,area_id)
            while !isempty(aq)
                root_area_id = pop_last!(aq)
                root_area = areas[root_area_id]
                
                if cond(root_area)
                    fun(root_area,area_id)
                    continue
                end

                for child_id in root_area.childs
                    push!(aq,child_id)
                end
            end
        end
    end
end





"""
    transform_topology_planar!(topo::Topology{3}, distort_fun, is_boundary)

Distort a structured hexahedral topology in-place while keeping all faces planar.

Applies separable, axis-aligned shifts:
    (x, y, z) -> (x + δx[x_layer], y + δy[y_layer], z + δz[z_layer])

Each δ is constant per coordinate layer; boundary nodes remain unchanged.

Arguments
- `topo::Topology{3}`: 3D topology (e.g., Cartesian hex grid)
- `distort_fun`: function called to produce an offset vector. Only components
  aligned with the respective axis are used to fill δx, δy, δz. The function
  is called with a single argument describing the layer: `(:x, xv)`, `(:y, yv)`,
  or `(:z, zv)`; implementations can ignore it.
- `is_boundary`: predicate `is_boundary(x)::Bool` to keep boundary nodes fixed.

Planarity rationale
- For an x-constant face, all nodes share the same `δx` (uniform translation),
  while `δy` and `δz` vary only along in-plane directions y and z.
- Analogous arguments hold for y- and z-constant faces.
"""
function transform_topology_planar!(
    topo::Topology{3},
    distort_fun::F1,
    is_boundary::F2,
) where {F1<:Function,F2<:Function}

    nodes = get_nodes(topo)
    n = length(nodes)
    n == 0 && return topo

    # Collect unique layer coordinates along each axis
    xs = unique([nodes[i][1] for i in 1:n]); sort!(xs)
    ys = unique([nodes[i][2] for i in 1:n]); sort!(ys)
    zs = unique([nodes[i][3] for i in 1:n]); sort!(zs)

    x_to_idx = Dict(x => i for (i, x) in enumerate(xs))
    y_to_idx = Dict(y => i for (i, y) in enumerate(ys))
    z_to_idx = Dict(z => i for (i, z) in enumerate(zs))

    δx = zeros(Float64, length(xs))
    δy = zeros(Float64, length(ys))
    δz = zeros(Float64, length(zs))

    # Fill axis-wise offsets from distort_fun; only use the aligned component
    for (i, xv) in enumerate(xs)
        off = distort_fun((:x, xv))
        v = SVector{3,Float64}(off)
        δx[i] = v[1]
    end
    for (j, yv) in enumerate(ys)
        off = distort_fun((:y, yv))
        v = SVector{3,Float64}(off)
        δy[j] = v[2]
    end
    for (k, zv) in enumerate(zs)
        off = distort_fun((:z, zv))
        v = SVector{3,Float64}(off)
        δz[k] = v[3]
    end

    # Apply to non-boundary nodes only
    for i in 1:n
        p = get_coords(nodes[i])
        is_boundary(p) && continue
        xi = x_to_idx[p[1]]
        yi = y_to_idx[p[2]]
        zi = z_to_idx[p[3]]
        newp = SVector{3,Float64}(p[1] + δx[xi], p[2] + δy[yi], p[3] + δz[zi])
        nodes[i] = Node(get_id(nodes[i]), newp)
    end

    return topo
end

"""
    transform_topology_linear_elements!(topo::Topology{3}, make_A, is_boundary)

Apply linear, per-element distortions while preserving face planarity and avoiding
conflicting node updates. Elements that share any node with an already-selected
one are skipped (node-disjoint selection). Elements that contain boundary nodes
(per `is_boundary`) are also skipped.

For each selected volume with node set V and barycenter c, the mapping is
x -> x + A*(x - c), where A = `make_A(c, diameter, vol_id)`.

- Faces remain planar because a linear map preserves planarity.
- Node-disjoint selection prevents two different maps from touching the same node.
- The diameter is computed as the maximum pairwise distance between nodes in V.

Arguments
- `topo::Topology{3}`: 3D topology to modify in-place
- `make_A`: function `(center::SVector{3,Float64}, diameter::Float64, vol_id::Int32) -> 3x3 matrix`
- `is_boundary`: predicate `is_boundary(x)::Bool` identifying boundary nodes
"""
function transform_topology_linear_elements!(
    topo::Topology{3},
    make_A::F1,
    is_boundary::F2,
) where {F1<:Function,F2<:Function}

    nodes = get_nodes(topo)
    n_nodes = length(nodes)
    n_nodes == 0 && return topo

    vol_node_ids_vec = get_volume_node_ids(topo)
    n_vols = length(vol_node_ids_vec)

    # Precompute candidate volumes: skip any volume containing a boundary node
    is_candidate = trues(n_vols)
    for vid in 1:n_vols
        for nid in vol_node_ids_vec[vid]
            if is_boundary(get_coords(nodes[nid]))
                is_candidate[vid] = false
                break
            end
        end
    end

    # Build node adjacency via edges: neighbors[i] lists nodes sharing an edge with i
    edge_nodes = get_edge_node_ids(topo)
    neighbors = [Int32[] for _ in 1:n_nodes]
    for en in edge_nodes
        a = en[1]; b = en[2]
        push!(neighbors[a], b)
        push!(neighbors[b], a)
    end

    # Select a subset with node and edge independence (no edge has both endpoints used)
    used = falses(n_nodes)
    selected = Int32[]
    for vid in 1:n_vols
        is_candidate[vid] || continue
        conflict = false
        for nid in vol_node_ids_vec[vid]
            if used[nid]
                conflict = true
                break
            end
            # Edge constraint: none of nid's neighbors are used
            for m in neighbors[nid]
                if used[m]
                    conflict = true
                    break
                end
            end
            conflict && break
        end
        if !conflict
            push!(selected, vid)
            for nid in vol_node_ids_vec[vid]
                used[nid] = true
            end
        end
    end

    # Helper to compute barycenter and diameter for a volume's nodes
    function _center_and_diameter(node_ids::AbstractVector{Int32})
        c = zero(SVector{3,Float64})
        for nid in node_ids
            c += get_coords(nodes[nid])
        end
        c /= length(node_ids)
        maxd2 = 0.0
        for i in 1:length(node_ids)
            pi = get_coords(nodes[node_ids[i]])
            for j in i+1:length(node_ids)
                pj = get_coords(nodes[node_ids[j]])
                d2 = sum(abs2, pi - pj)
                if d2 > maxd2
                    maxd2 = d2
                end
            end
        end
        return c, sqrt(maxd2)
    end

    # Apply linear map per selected volume
    for vid in selected
        node_ids = vol_node_ids_vec[vid]
        c, diam = _center_and_diameter(node_ids)
        A = make_A(c, diam, vid)
        @assert size(A) == (3,3) "make_A must return a 3x3 matrix"
        for nid in node_ids
            p = get_coords(nodes[nid])
            # Safety: boundary nodes should not be present, but double-check
            is_boundary(p) && continue
            newp = p + A*(p - c)
            nodes[nid] = Node(nid, newp)
        end
    end

    return topo
end

"""
    transform_topology_linear_elements!(topo::Topology{3}, is_boundary; fraction=1/8)

Default variant that uses an isotropic scaling A = s*I per selected element so
that farthest nodes move by approximately `fraction*diameter` (for cubes, s=2*fraction).
"""
function transform_topology_linear_elements!(
    topo::Topology{3},
    is_boundary::F2;
    fraction::Float64 = 1/8,
) where {F2<:Function}
    s = 2*fraction
    make_A = (_c, _diam, _vid) -> SMatrix{3,3,Float64}(I) * s
    return transform_topology_linear_elements!(topo, make_A, is_boundary)
end
