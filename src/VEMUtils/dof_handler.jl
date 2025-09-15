
struct DofHandler{D,U}
    dof_mapping::Dict{Int,SVector{U,Int}}
end


function DofHandler{U}(mesh::Mesh{D}) where {D,U}

    counter = 1
    dof_mapping = Dict{Int,SVector{U,Int}}()
    for node in mesh.nodes
        node_id = get_id(node)
        !is_active(node) && continue
        dof_mapping[node_id] = SVector{U,Int}(counter:counter+U-1)
        counter += U
    end

    return DofHandler{D,U}(dof_mapping)
end



@inline get_dofs(dh::DofHandler,id::Integer) = dh.dof_mapping[id]

"""
    get_dof(dh::DofHandler{U},id::Int,i::Int) where {U}

id is the node_id and i is the i-th component of the vector field
"""
function get_dof(dh::DofHandler,id::Integer,i::Integer)
    return get_dofs(dh,id)[i]
end


function get_dofs!(
    vec::AbstractVector{Int},
    dh::DofHandler{D,U},
    nodes::AbstractVector{N}; 
    start = 1) where {D,U,N}

    n_dofs = length(nodes)*U + (start-1)
    length(vec) != n_dofs && resize!(vec,n_dofs)

    counter = start
    @inbounds for node in nodes,i in 1:U
        id = get_id(node)  # works for nodes and integers
        vec[counter] = get_dof(dh,id,i)
        counter += 1
    end
end

function get_dofs(dh::DofHandler{U},
    nodes::AbstractVector) where U

    vec = Vector{Int}(undef,length(nodes)*U)
    get_dofs!(vec,dh,nodes)
    return vec
end