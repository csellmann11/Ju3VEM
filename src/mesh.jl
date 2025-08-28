using FastGaussQuadrature
using Ju3VEM
using Ju3VEM.VEMGeo: get_base_len, FaceIntegralData,resize_and_fill!,get_outward_normal
using Ju3VEM.VEMGeo: get_iterative_area_vertex_ids!,get_edge,get_area,get_bc,get_hf
using Ju3VEM.VEMGeo: get_gauss_lobatto,get_normal,get_exp_to_idx_dict,precompute_volume_monomials
using FixedSizeArrays,LinearAlgebra
using StaticArrays,LoopVectorization
using Octavian
using Chairmarks

abstract type ElType{K} end

struct StandardEl{K} <: ElType{K} end
struct SerendipityEl{K} <: ElType{K} end


#TODO: remove this, already in VEMGeo
const _GaussLobattoLookup = [gausslobatto(i) for i in 2:8]
gauss_lobatto(i::Int) = _GaussLobattoLookup[i]

@kwdef struct Mesh{D,ET}
    nodes::FlattenVecs{4,Node{D,Float64},Vector{Node{D,Float64}}}

    topo::Topology{D}
    int_coords_connect::Vector{Vector{Vector{Int}}} #TODO: Improve naming

    node_sets::Dict{String,Set{Int}} = Dict{String,Set{Int}}()
    edge_sets::Dict{String,Set{Int}} = Dict{String,Set{Int}}()
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

    #INFO: For VEM and HHO, g_nodes are used for integration, therfore they cant be equally spaced
    g_nodes = gauss_lobatto(K)[1][2:end-1]

    edges = topo.connectivity[1, 2]
    for edge in RootIterator{D,2}(topo)
        c1, c2 = vertices[edges[edge.id]]
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
        node_ids = get_volume_node_ids(topo, vol.id)
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



###############################################################
# Move the following to a test file
###############################################################

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




mesh = Mesh(topo, StandardEl{1}())

n_edges = length(get_edges(topo))

get_vertices(mesh)
get_gauss_nodes(mesh)
get_face_moments(mesh)
get_volume_moments(mesh)



include("face_projector.jl")
include("element_mapping.jl")
include("volume_projector.jl")

base   = get_base(BaseInfo{2,1,1}())

function create_facedata_col(mesh::Mesh{3,StandardEl{K}}) where K

    K_MAX = max(2*K-1,2)
    base = get_base(BaseInfo{2,K_MAX,1}()).base
    topo = mesh.topo
    facedata_col = Dict{Int,FaceData{3,length(base)}}()
    # volume integrals are needed up to order max(2*(K-1),1)
    # -> since the divergence theorem is used we need the 
    # face integrals up to order max(2*(K-1),1) + 1 = max(2*K-1,2)


    for face in RootIterator{3}(topo)
        dΩ = precompute_face_monomials(face.id, topo, Val(K_MAX))
        facedata_col[face.id] = h1_projectors!(face.id,mesh,dΩ)
    end
    return facedata_col
end

facedata_col = create_facedata_col(mesh)

ntl = create_node_mapping(1,mesh,facedata_col)
vol_id = 1


vol_data = precompute_volume_monomials(1,topo,facedata_col,Val(1))

Ju3VEM.VEMGeo.compute_volume_integral_unshifted(Monomial(1.0,SA[0,0,0]),vol_data)
bmat3d = create_volume_bmat(vol_id,mesh,vol_data.vol_bc,vol_data.hvol,facedata_col,ntl)

@b create_volume_bmat($vol_id,$mesh,$vol_data.vol_bc,$vol_data.hvol,$facedata_col,$ntl)






# A = rand(4,5); 
# B = rand(5,6); 

# C1 = static_matmul(A,B,Val((4,6)))
# C2 = A*B
# C3 = static_matmul_turbo(A,B,Val((4,6)))

# C1 ≈ C2 ≈ C3

# @b static_matmul($A,$B,Val((4,6)))
# @b static_matmul_turbo($A,$B,Val((4,6)))
# @b $A*$B