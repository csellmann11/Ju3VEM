using FixedSizeArrays,LinearAlgebra
using StaticArrays,LoopVectorization
using ..VEMGeo: get_base_len, FaceIntegralData, resize_and_fill!, get_outward_normal
using ..VEMGeo: get_iterative_area_vertex_ids!, get_edge, get_area, get_bc, get_hf
using ..VEMGeo: get_gauss_lobatto, get_normal, get_exp_to_idx_dict, precompute_volume_monomials
using ..VEMGeo: interpolate_edge

abstract type ElType{K} end

struct StandardEl{K} <: ElType{K} end
struct SerendipityEl{K} <: ElType{K} end


# Use VEMGeo.get_gauss_lobatto

@kwdef struct Mesh{D,ET}
    nodes::FlattenVecs{4,Node{D,Float64},Vector{Node{D,Float64}}}

    topo::Topology{D}
    int_coords_connect::Vector{Vector{Vector{Int}}} #TODO: Improve naming

    node_sets::Dict{String,Set{Int}} = Dict{String,Set{Int}}()
    edge_sets::Dict{String,Set{Int}} = Dict{String,Set{Int}}()
    face_sets::Dict{String,Set{Int}} = Dict{String,Set{Int}}()
    volume_sets::Dict{String,Set{Int}} = Dict{String,Set{Int}}()
end


function Mesh(topo::Topology{D}, ::ET) where {D,ET<:ElType}

    # finish_topology!(topo) #TODO: write this code

    node_type = eltype(get_nodes(topo))

    int_coords_connect = Vector{Vector{Vector{Int}}}(undef, D + 1)

    int_coords_connect[2] = [Int[] for _ in get_edges(topo)]
    int_coords_connect[3] = [Int[] for _ in get_areas(topo)]
    if D >= 3
        int_coords_connect[4] = [Int[] for _ in get_volumes(topo)]
    end

    nodes = FlattenVecs(get_nodes(topo), node_type[], node_type[], node_type[])

    mesh = Mesh{D,ET}(nodes=nodes, topo=topo,
        int_coords_connect=int_coords_connect)

    add_internal_coords!(mesh, int_coords_connect)

    mesh
end


@inline get_order(::Mesh{D,ET}) where {D,K,ET<:ElType{K}} = K

Base.@propagate_inbounds Base.getindex(mesh::Mesh,idx::Integer) = mesh.nodes[idx]


get_vertices(mesh::Mesh) = mesh.nodes.v.args[1]
get_gauss_nodes(mesh::Mesh) = mesh.nodes.v.args[2]
get_face_moments(mesh::Mesh) = mesh.nodes.v.args[3]
get_volume_moments(mesh::Mesh) = mesh.nodes.v.args[4]


function getnfacemoments(::Int, ::StandardEl{K}, topo::Topology{D}) where {D,K}
    base_len = get_base_len(2, K - 2, 1)
    return base_len
end

function getnvolmoments(::Int, ::StandardEl{K}, topo::Topology{D}) where {D,K}
    base_len = get_base_len(3, K - 2, 1)
    return base_len
end

function getnfacemoments(area_id::Int, 
    ::SerendipityEl{K}, topo::Topology{D}; lazy=true) where {D,K}
    if lazy
        return get_base_len(2, K - 3, 1)
    end
    geo_edge_ids = get_area_edge_ids(topo, area_id)
    ηE = length(geo_edge_ids)
    # base_order = K - ηE + 2 #TODO: why did we add 2?
    base_order = K - ηE
    return get_base_len(2, base_order, 1)
end

function getnvolmoments(vol_id::Int, 
    ::SerendipityEl{K}, topo::Topology{D}; lazy=true) where {D,K}
    if lazy
        return get_base_len(3, K - 3, 1)
    end
    geo_area_ids = get_volume_area_ids(topo, vol_id)
    ηA = length(geo_area_ids)
    # base_order = K - ηA + 2 #TODO: why did we add 2?
    base_order = K - ηA
    return get_base_len(3, base_order, 1)
end


function get_iterative_area_node_ids!(
    face_node_ids::FlattenVecs{3,Int},
    area::Area{D}, mesh::Mesh{D,ET},
    abs_ref_level::Int=typemax(Int)) where {D,K,ET<:ElType{K}}

    topo = mesh.topo
    vertex_ids = face_node_ids.v.args[1]
    get_iterative_area_vertex_ids!(vertex_ids, area, topo, abs_ref_level)
    K == 1 && return nothing

    n_edge_vertices = length(vertex_ids)*(K-1) 
    edge_vertices   = face_node_ids.v.args[2]
    resize!(edge_vertices,n_edge_vertices)

    for i in eachindex(vertex_ids)
        ip1 = get_next_idx(vertex_ids,i)
        vertex1 = vertex_ids[i]
        vertex2 = vertex_ids[ip1]
        edge = get_edge(vertex1,vertex2,topo)
        
        for (j,v) in enumerate(mesh.int_coords_connect[2][edge.id])
            edge_vertices[(i-1)*(K-1) + j] = v
        end
    end

    moment_ids = face_node_ids.v.args[3]
    resize_and_fill!(moment_ids,mesh.int_coords_connect[3][area.id])
    return nothing
end

function add_internal_coords!(mesh::Mesh{D,ET},
    int_coords_connect::Vector{Vector{Vector{Int}}}
            ) where {D,K,ET<:ElType{K}}

    topo = mesh.topo

    gauss_coords = get_gauss_nodes(mesh)
    face_moment_coords = get_face_moments(mesh)
    volume_moment_coords = get_volume_moments(mesh)
    empty!(gauss_coords)
    empty!(face_moment_coords)
    empty!(volume_moment_coords)
    # if O == 1 return nothing end


    vertices = get_nodes(topo)
    num_mesh_nodes = length(get_vertices(mesh))

    #INFO: For VEM and HHO, g_nodes are used for integration, therefore they can't be equally spaced
    # K is the polynomial order; use K+1 Lobatto nodes and drop endpoints → interior edge nodes
    g_nodes = get_gauss_lobatto(K+1)[1][2:end-1]

    # edges = topo.connectivity[1, 2]
    edge_nodes = get_edge_node_ids(topo)
    for edge in RootIterator{D,2}(topo)
        c1, c2 = vertices[edge_nodes[edge.id]]
        # edge_gauss_ids = int_coords_connect[2][edge.id]

        edge_gauss_ids = Vector{Int}(undef, length(g_nodes))
        for i in eachindex(g_nodes)
            new_point_coords = interpolate_edge(g_nodes[i], c1, c2)

            new_id = length(gauss_coords) + num_mesh_nodes + 1
            new_node = Node(new_id, new_point_coords)

            #TODO: add boundary info
            # push!(topo.boundary_vec,
            #     is_boundary_edge(topo, edge))


            push!(gauss_coords, new_node) # push node
            edge_gauss_ids[i] = new_id
        end
        int_coords_connect[2][edge.id] = edge_gauss_ids
    end

    num_edge_vertices = length(gauss_coords)
    face_moment_offset = num_edge_vertices + num_mesh_nodes

    for area in RootIterator{D,3}(topo)
        #Info: Moments get a coord for later dof reordering
        # node_ids = topo.connectivity[1, 3][area.id]
        node_ids = get_area_node_ids(topo, area.id)

        num_moments_local = getnfacemoments(area.id, ET(), topo)
        area_moment_ids = FixedSizeVector{Int}(undef, num_moments_local)

        for i in eachindex(area_moment_ids)
            # new_point_coords = mean(vertices[id] for id in node_ids)
            new_point_coords = zero(Node{D,Float64})
            new_id = length(face_moment_coords) + face_moment_offset + 1
            new_node = Node(new_id, new_point_coords)
            push!(face_moment_coords, new_node)
            area_moment_ids[i] = new_id
        end

        int_coords_connect[3][area.id] = area_moment_ids
    end


    volume_moment_offset = face_moment_offset + length(face_moment_coords)
    for vol in RootIterator{D,4}(topo)
        # node_ids = get_volume_node_ids(topo, vol.id) 
        num_moments_local = getnvolmoments(vol.id, ET(), topo)
        volume_moment_ids = FixedSizeVector{Int}(undef, num_moments_local)

        for i in eachindex(volume_moment_ids)
            new_point_coords = zero(Node{D,Float64})
            new_id = length(volume_moment_coords) + volume_moment_offset + 1
            new_node = Node(new_id, new_point_coords)
            push!(volume_moment_coords, new_node)
            volume_moment_ids[i] = new_id
        end
        int_coords_connect[4][vol.id] = volume_moment_ids
    end
end



function add_mesh_nodes_unique(topo, 
    nodes_added::Set{Int})

    pot_id = length(topo.nodes) + 1
    if pot_id ∉ nodes_added
        push!(nodes_added, pot_id)
        return pot_id
    end
    return -1
end




