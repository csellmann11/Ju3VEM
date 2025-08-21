using FastGaussQuadrature


abstract type ElType{K} end

struct StandardEl{K} <: ElType{K} end 
struct SerendipityEl{K} <: ElType{K} end


const _GaussLobattoLookup = [gausslobatto(i) for i in 2:8]
gauss_lobatto(i::Int) = _GaussLobattoLookup[i]

@kwdef struct Mesh{D,ET}
    nodes  ::FlattenVecs{D+1,Node{D},Vector{Node{D}}}

    topo              ::Topology{D}
    int_coords_connect::Vector{Vector{Vector{Int}}} #TODO: Improve naming

    node_sets   ::Dict{String,Set{Int}} = Dict{String,Set{Int}}()
    edge_sets   ::Dict{String,Set{Int}} = Dict{String,Set{Int}}()
end


function Mesh(topo::Topology{D},::ET) where {D,ET<:ElType}
    
    # finish_topology!(topo) #TODO: write this code

    node_type = eltype(get_nodes(topo))

    int_coords_connect = Vector{Vector{Vector{Int}}}(undef,D+1)

    int_coords_connect[2] = [Int[] for _ in get_edges(topo)]
    int_coords_connect[3] = [Int[] for _ in get_areas(topo)]
    if D >= 3 
        int_coords_connect[4] = [Int[] for _ in get_volumes(topo)]
    end

    nodes = FlattenVecs(get_nodes(topo),node_type[],node_type[])

    mesh =  Mesh{D,ET}(nodes = nodes, topo = topo,
        int_coords_connect = int_coords_connect)

    add_internal_coords!(mesh,int_coords_connect) 
    
    mesh
end


get_vertices(mesh::Mesh)       = mesh.nodes.v.args[1]
get_gauss_nodes(mesh::Mesh)    = mesh.nodes.v.args[2]
get_face_moments(mesh::Mesh)   = mesh.nodes.v.args[3]
get_volume_moments(mesh::Mesh) = mesh.nodes.v.args[4]


function getnfacemoments(::Int,::StandardEl{K},topo::Topology{D}) where {D,K}
    base_len = get_base_len(D,K-2,1)
    return base_len
end

function getnfacemoments(area_id::Int,::SerendipityEl{K},topo::Topology{D};lazy=true) where {D,K}
    if lazy 
        return get_base_len(D,K-3,1)
    end
    geo_edge_ids = get_area_edge_ids(topo,area_id)
    ηE = length(geo_edge_ids) 
    base_order = K-ηE+2
    return get_base_len(D,base_order,1)
end


function add_internal_coords!(mesh::Mesh{D,ET},
    int_coords_connect::Vector{Vector{Vector{Int}}}
    ) where {D,K,ET<:ElType{K}}

topo = mesh.topo

gauss_coords = gauss_nodes(mesh)
moment_coords = moments(mesh)
empty!(gauss_coords); empty!(moment_coords)
# if O == 1 return nothing end


vertices = get_nodes(topo)
num_mesh_nodes  = length(mesh_nodes(mesh))

#INFO: For VEM and HHO, g_nodes are used for integration, therfore they cant be equally spaced
g_nodes = gauss_lobatto(K+1)[1][2:end-1]

edges = topo.connectivity[1,2]
for edge in RootIterator{D,2}(topo)
    c1,c2 = vertices[edges[edge.id]]
    # edge_gauss_ids = int_coords_connect[2][edge.id]

    edge_gauss_ids = Vector{Int}(undef,length(g_nodes))
    for i in eachindex(g_nodes)
        new_point_coords = interpolate_edge(g_nodes[i],c1,c2)

        new_id = length(gauss_coords) + num_mesh_nodes + 1
        new_node = Node{D}(new_id,new_point_coords)

        push!(topo.boundary_vec,
                    is_boundary_edge(topo,edge)) 
   

        push!(gauss_coords,new_node) # push node
        edge_gauss_ids[i] = new_id
    end
    int_coords_connect[2][edge.id] = edge_gauss_ids
end

num_gauss_nodes = length(gauss_coords)

for area in RootIterator{D,3}(topo)
    #Info: Moments get a coord for later dof reordering
    node_ids        = topo.connectivity[1,3][area.id]

    num_moments_local = get_num_moments(ElT,topo,area.id) 
    area_moment_ids = Vector{Int}(undef,num_moments_local)

    for i in eachindex(area_moment_ids)
        new_point_coords = mean(vertices[id] for id in node_ids) 
        new_id = length(moment_coords) + num_mesh_nodes + num_gauss_nodes + 1
        new_node = Node{D}(new_id,new_point_coords)
        push!(moment_coords,new_node)
        area_moment_ids[i] = new_id
    end

    int_coords_connect[3][area.id] = area_moment_ids
end
end