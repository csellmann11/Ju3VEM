

mutable struct Node{D,T} <: StaticVector{D,T}
    id::Int
    const coords::SVector{D,T}
end

function Node(id::Int, coords::AbstractVector{T}) where T
    D = length(typeof(coords))
    Node{D,T}(id, coords)
end

Node(coords::AbstractVector) = Node(-1, coords)
Node(coords::StaticVector) = Node(-1, coords)



get_coords(n::Node) = n.coords
get_id(n::Node) = n.id
get_id(i::Integer) = i
Base.getindex(n::Node{D,T}, i::Int) where {D,T} = get_coords(n)[i]

Base.zero(::Type{Node{D,T}}) where {D,T} = Node(zero(SVector{D,T})) # to create placeholder nodes
Base.zero(n::Node) = zero(typeof(n))
Base.rand(::Type{Node{D,T}}) where {D,T} = Node(-1, @SVector rand(T,D))


#########################
# NManifold
#########################
@kwdef mutable struct NManifold{N,D}
    id::Int
    is_root::Bool = true
    const refinement_level::Int

    # tree structure
    const parent_id::Int
    const childs::Vector{Int} = Int[]
end

const Edge{D} = NManifold{2,D}
const Area{D} = NManifold{3,D}
const Volume{D} = NManifold{4,D}

function NManifold{N,D}(id::Int, refinement_level::Int, parent_id::Int) where {N,D}
    NManifold{N,D}(id=id, refinement_level=refinement_level, parent_id=parent_id)
end

function activate!(mf::MF) where MF <: Union{<:NManifold,<:Node}
    mf.id = abs(mf.id)
    return mf
end

function deactivate!(mf::Union{<:NManifold,<:Node})
    mf.id = -abs(mf.id)
    return mf
end

get_id(mf::NManifold) = mf.id
is_active(m::Union{Node,NManifold}) = m.id > 0


#########################
# Topology
#########################
using FixedSizeArrays
Base.convert(::FixedSizeVector{T,V},v::SVector{N,T}) where {T,N,V} = FixedSizeVector{T,V}(v) 
const ID_VEC_TYPE{T} = FixedSizeVector{T,Memory{T}}
@kwdef struct Topology{D}

    nodes::Vector{Node{D,Float64}}              = Node{D,Float64}[]
    edges::OrderedDict{ID_VEC_TYPE{Int},Edge{D}}     = OrderedDict{ID_VEC_TYPE{Int},Edge{D}}()
    areas::OrderedDict{ID_VEC_TYPE{Int},Area{D}}     = OrderedDict{ID_VEC_TYPE{Int},Area{D}}()
    volumes::OrderedDict{ID_VEC_TYPE{Int},Volume{D}} = OrderedDict{ID_VEC_TYPE{Int},Volume{D}}()


    connectivity::Matrix{Vector{ID_VEC_TYPE{Int}}} = begin
        mat = Matrix{Vector{ID_VEC_TYPE{Int}}}(undef, 4, 4)
        mat[1, 2] = edges.keys
        mat[1, 3] = areas.keys
        mat[3, 4] = volumes.keys #! volumes are index by face_ids
        mat[1, 4] = Vector{ID_VEC_TYPE{Int}}()
        for i in 1:D+1, j in 2:i-1
            if i != j && (i,j) != (4,3)
                mat[j, i] = Vector{ID_VEC_TYPE{Int}}()
            end
        end 
        mat
    end

    el_neighs::Dict{Int,Vector{Int}} = Dict{Int,Vector{Int}}()

    #Here also hanging nodes are contained wich are not found in the connectivity
    nids_col::Dict{Int,Vector{Int}} = Dict{Int,Vector{Int}}()
    boundary_vec::BitVector = falses(0)
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
@inline get_volume_node_ids(topo)::Vector{ID_VEC_TYPE{Int}}      = topo.connectivity[1, 4]
@inline get_volume_node_ids(topo, id::Int)::ID_VEC_TYPE{Int}     = get_volume_node_ids(topo)[id]
@inline get_volume_edge_ids(topo)::Vector{ID_VEC_TYPE{Int}}      = topo.connectivity[2, 4]
@inline get_volume_edge_ids(topo, id::Int)::ID_VEC_TYPE{Int}     = get_volume_edge_ids(topo)[id]
@inline get_volume_area_ids(topo)::Vector{ID_VEC_TYPE{Int}}      = topo.connectivity[3, 4]
@inline get_volume_area_ids(topo, id::Int)::ID_VEC_TYPE{Int}     = get_volume_area_ids(topo)[id]
#Info: access face information
@inline get_area_node_ids(topo)::Vector{ID_VEC_TYPE{Int}}        = topo.connectivity[1, 3]
@inline get_area_node_ids(topo, id::Int)::ID_VEC_TYPE{Int}       = get_area_node_ids(topo)[id]
@inline get_area_edge_ids(topo)::Vector{ID_VEC_TYPE{Int}}        = topo.connectivity[2, 3]
@inline get_area_edge_ids(topo, id::Int)::ID_VEC_TYPE{Int}       = get_area_edge_ids(topo)[id]
#Info: access edge information
@inline get_edge_node_ids(topo)::Vector{ID_VEC_TYPE{Int}}        = topo.connectivity[1, 2]
@inline get_edge_node_ids(topo, id::Int)::ID_VEC_TYPE{Int}       = get_edge_node_ids(topo)[id]


@inline add_volume_node_id!(topo, dest::Int, val::Int) =
    push!(get_volume_node_ids(topo, dest), val)
@inline add_volume_edge_id!(topo, dest::Int, val::Int) =
    push!(get_volume_edge_ids(topo, dest), val)
@inline add_volume_area_id!(topo, dest::Int, val::Int) =
    push!(get_volume_area_ids(topo, dest), val)

@inline add_area_node_id!(topo, dest::Int, val::Int) =
    push!(get_area_node_ids(topo, dest), val)
@inline add_area_node_id!(topo,loc::Int,dest::Int,val::Int) = 
    get_area_node_ids(topo,dest)[loc] = val


@inline add_area_edge_id!(topo, dest::Int, val::Int) =
    push!(get_area_edge_ids(topo, dest), val)
@inline add_area_edge_id!(topo,loc::Int,dest::Int,val::Int) = 
    get_area_edge_ids(topo,dest)[loc] = val


@inline add_edge_node_id!(topo, dest::Int, val::Int) =
    push!(get_edge_node_ids(topo, dest), val)
@inline add_edge_node_id!(topo,loc::Int,dest::Int,val::Int) = 
    get_edge_node_ids(topo,dest)[loc] = val


@inline __default_node_id(topo) = length(get_nodes(topo)) + 1


function add_node!(coords::AbstractVector,
    topo::Topology{D},
    id=__default_node_id(topo)) where {D}

    node = Node(id, coords)
    push!(get_nodes(topo), node)
    return id
end
function add_node!(node::Node,
    topo::Topology)
    add_node!(node.coords, topo, node.id)

    return node.id
end


function rotate_to_min(v::AbstractVector{<:Integer})
    # Rotate the vector so that the smallest element is first
    dest = ID_VEC_TYPE{Int}(undef, length(v))
    circshift!(dest,v,-(argmin(v) - 1))
    return dest
end

@inline edge_hash(n1::Int, n2::Int) = ID_VEC_TYPE(SVector(minmax(n1, n2)))
@inline edge_hash(v::Union{SVector{2,Int},Tuple{Int,Int}}) = edge_hash(v...)
@inline area_hash(v) = rotate_to_min(v)
@inline volume_hash(v) = sort!(ID_VEC_TYPE(v))

function add_edge!(node_ids::AbstractVector,
    topo::Topology{D},
    parent_id::Int=0,
    parent_ref_level::Int=0) where D

    @assert length(node_ids) == 2

    temp_id = length(topo.edges) + 1
    edge = get!(topo.edges, edge_hash(node_ids),
        Edge{D}(temp_id, parent_ref_level + 1, parent_id))

    activate!(edge)
    return edge.id
end

function activate!(area::Area{D}, topo::Topology{D})::Nothing where D
    area.id > 0 && return nothing
    activate!(area)
    for edge_id in get_area_edge_ids(topo,area.id)
        edge = get_edges(topo)[edge_id]
        activate!(edge)
    end
    nothing
end

function add_area!(_node_ids::AbstractVector{Int},
    topo::Topology{D},
    parent_id::Int=0,
    parent_ref_level::Int=0) where D

    area_dict = topo.areas
    temp_id = length(area_dict) + 1


    area = get!(area_dict, area_hash(_node_ids),
        Area{D}(temp_id, parent_ref_level + 1, parent_id))

 
    !is_active(area) && activate!(area, topo)
    id = area.id
    id != temp_id && return area.id

    node_ids = get_area_node_ids(topo, id)


    area_edge_ids = ID_VEC_TYPE{Int}(undef,length(node_ids))
    push!(get_area_edge_ids(topo), area_edge_ids)
    for (i, node_id) in enumerate(node_ids)
        ip1 = get_next_idx(node_ids, i)
        edge_id = add_edge!(SA[node_id, node_ids[ip1]],
            topo, 0, parent_ref_level)
 
        add_area_edge_id!(topo,i,id,edge_id)
    end
    return id
end

function activate!(volume::Volume{D}, topo::Topology{D})::Nothing where D
    volume.id > 0 && return nothing
    activate!(volume)

    for area_id in get_volume_area_ids(topo,volume.id)
        area = get_areas(topo)[area_id]
        activate!(area)
    end
    nothing
end


function add_volume!(face_ids::AbstractVector{Int},
    topo::Topology{D},
    parent_id::Int=0,
    parent_ref_level::Int=0) where D
    
    volume_dict = topo.volumes
    temp_id = length(volume_dict) + 1

    volume = get!(volume_dict, volume_hash(face_ids),
        Volume{D}(temp_id, parent_ref_level + 1, parent_id))

    !is_active(volume) && activate!(volume, topo)
    id = volume.id 

    
    if id == temp_id
        unique_ids = get_unique_values(get_area_node_ids(topo),face_ids)
        push!(get_volume_node_ids(topo), unique_ids)
    end
    volume.id

end


function add_volume!(node_ids_col::AbstractVector{<:AbstractVector{Int}},
    topo::Topology{D},
    parent_id::Int=0,
    parent_ref_level::Int=0) where D

    n_faces = length(node_ids_col)
    volume_dict = topo.volumes
    temp_id = length(volume_dict) + 1
 

    _face_ids = ID_VEC_TYPE{Int}(undef, n_faces)

    for (i, node_ids) in enumerate(node_ids_col)
        face_id = add_area!(node_ids, topo, 0, parent_ref_level)
        _face_ids[i] = face_id
    end

    volume = get!(volume_dict, volume_hash(_face_ids),
        Volume{D}(temp_id, parent_ref_level + 1, parent_id))

    activate!(volume, topo)
    id = volume.id


    if id == temp_id
        push!(get_volume_node_ids(topo), get_unique_values(node_ids_col))
    end
    volume.id
    
end




@inline is_root(mf::NManifold) = mf.is_root
@inline is_active_root(mf::NManifold) = is_root(mf) && is_active(mf)

#TODO: check performance implications of union return type
@inline function get_edge(n1::Int, n2::Int, topo::Topology{D}) where D
    edge = get(topo.edges, edge_hash(SVector(n1, n2)), nothing)
    return edge
end

@inline function has_edge(n1::Int, n2::Int, topo::Topology{D}) where D
    edge_hash = edge_hash(SVector(n1, n2))
    return haskey(topo.edges, edge_hash)
end
@inline has_edge(n1::Node, n2::Node, topo::Topology) = has_edge(n1.id, n2.id, topo)

@inline get_edge(n1::Node, n2::Node, topo::Topology) = get_edge(n1.id, n2.id, topo)

@inline function get_face(v::AbstractVector{Int}, topo::Topology{D}) where D
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
    state_is_root::BitVector = is_root.(get_geo_vec(topo, Val(Idx)))
    return RootIterator{D,Idx}(topo, state_is_root)
end

function RootIterator{D,1}(topo::Topology{D}) where {D}
    # state_is_root::BitVector = is_root.(get_geo_vec(topo, Val(Idx)))
    state_is_root = trues(length(get_nodes(topo)))
    return RootIterator{D,1}(topo, state_is_root)
end

Base.firstindex(r::RootIterator{D,Idx}) where {D,Idx} = 1

function Base.iterate(r::RootIterator{D,Idx}, state=1) where {D,Idx}

    while state <= length(r.state_is_root)
        geo = get_geo_vec(r.topo, Val(Idx))[state]

        # is_active check is needed if an element has been coarsed 

        # When an element is coarsed also the other childs are coarsed
        # which are yet to be iterated
        if r.state_is_root[state] && is_active(geo)
            return geo, state + 1
        end
        state += 1
    end
    return nothing
end


function num_roots(geo_vec::AbstractVector{T}) where T<:NManifold
    return count(geo -> is_root(geo) && is_active(geo), geo_vec)
end

# Base.length(r::RootIterator{D,Idx}) where {D,Idx} = sum(r.state_is_root)
Base.length(r::RootIterator{D,Idx}) where {D,Idx} = num_roots(get_geo_vec(r.topo, Val(Idx)))


#############################################
# Tree utils
#############################################
"""
    get_iterative_area_vertex_ids(area::Area{D}, topo::Topology{D},
        abs_ref_level::Int=typemax(Int)) where D

Get the node ids of the mesh of an area.

# Arguments
- `area::Area{D}`: The area to get the mesh node ids of.
- `topo::Topology{D}`: The topology to get the mesh node ids from.
- `abs_ref_level::Int`: Does not go below this refinement level in the recursive search
"""
function get_iterative_area_vertex_ids!(
    vertex_ids::AbstractVector{Int},
    area::Area{D}, topo::Topology{D},
    abs_ref_level::Int=typemax(Int)) where D

    empty!(vertex_ids)
    cond(edge) = edge.refinement_level == abs_ref_level || is_root(edge)
    iterate_element_edges(topo, area.id, cond) do node_id, _, _
        push!(vertex_ids, node_id)
    end
    return nothing
end


function get_iterative_area_vertex_ids(area::Area{D}, topo::Topology{D},
    abs_ref_level::Int=typemax(Int)) where D

    vertex_ids = Int[]
    get_iterative_area_vertex_ids!(vertex_ids, area, topo, abs_ref_level)
    return vertex_ids
end

"""
iterate_element_edges(f::F1,topo::Topology{D},area_id::Int, cond::F2 = is_root)

Iterate over the edges of an element and apply a function on the edges

# Arguments
- `f::Function`: The function to apply on the edges. f takes 
the `node_id`, `edge_id` and `parent_edge_id` as arguments
"""
function iterate_element_edges(fun::F1, topo::Topology{D}, area_id::Int, cond::F2=is_root) where {D,F1,F2}
    edge_ids = topo.connectivity[2, 3][area_id]
    edges = get_edges(topo)

    nodes_of_edges = get_edge_node_ids(topo)
    area_nodes = get_area_node_ids(topo, area_id)

    for (count, parent_edge_id) in enumerate(edge_ids)
        edge = edges[parent_edge_id]

        n1 = area_nodes[count]
        countp1 = get_next_idx(area_nodes, count)
        n2 = area_nodes[countp1]

        apply_f_on(cond, edge, edges, n2 < n1) do root_edge, rev_child_order
            # idx = rev_child_order ? 1 : 2
            idx = 1 + rev_child_order
            edge_id = root_edge.id
            node_id = nodes_of_edges[edge_id][idx]
            fun(node_id, edge_id, parent_edge_id)
        end
    end
end


function iterate_volume_areas(fun::F1, topo::Topology{D}, volume_id::Int, cond::F2=is_root) where {D,F1,F2}
    area_ids = get_volume_area_ids(topo, volume_id)
    areas = get_areas(topo)
    for area_id in area_ids
        area = areas[area_id]
        apply_f_on(cond, area, areas, false) do root_area, _
            fun(root_area, area, volume_id, topo)
        end
    end
end




"""
    apply_f_on(f::F1, cond::F2,
        feature::T, features::AbstractVector{T},
        rev_child_order::Bool) where {T,F1<:Function,F2<:Function}

Apply a function on a feature and its children if the condition is true.

# Arguments 
- `f::Function`: The function to apply on the feature. f takes 
the `feature`, `rev_child_order` as arguments
- `cond::Function`: The condition to check if the feature should be applied.
- `feature::T`: The feature to apply the function on.
- `features::AbstractVector{T}`: The vector of features to apply the function on.
- `rev_child_order::Bool`: If true, the child order is reversed.
"""
function apply_f_on(f::F1, cond::F2,
    feature::T, features::AbstractVector{T},
    rev_child_order::Bool) where {T,F1<:Function,F2<:Function}
    if cond(feature) && is_active(feature)
        return f(feature, rev_child_order)
    end

    # If this feature has no children, terminate traversal safely
    if isempty(feature.childs)
        return nothing
    end

    child_id1 = feature.childs[1+rev_child_order]
    apply_f_on(f, cond, features[child_id1], features, false)

    child_id2 = feature.childs[2-rev_child_order]
    apply_f_on(f, cond, features[child_id2], features, true)
end

apply_f_on_roots(f, feature, features, rev_child_order::Bool=false) =
    apply_f_on(f, is_root, feature, features, rev_child_order)






