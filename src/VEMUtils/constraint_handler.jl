function check_unique_set_name(mesh::Mesh,set_name::String,set_type::Symbol)

    if set_type == :node_sets
        set = mesh.node_sets
    elseif set_type == :edge_sets
        set = mesh.edge_sets
    else
        error("Invalid set type")
    end

    if haskey(set,set_name)
        error("Set name $set_name already exists")
    end

    return true
end

function add_node_set!(mesh::Mesh,
    set_name::String,
    f::F;
    check_unique::Bool = true) where F <: Function

    check_unique && check_unique_set_name(mesh,set_name,:node_sets)

    n_set = Set{Int}()
    for node in Iterators.flatten((get_vertices(mesh),get_gauss_nodes(mesh)))

        !is_active(node) && continue
        f(node) && push!(n_set,get_id(node))
    end

    mesh.node_sets[set_name] = n_set
end


function add_edge_set!(mesh::Mesh,
    set_name::String,
    f::F;
    check_unique::Bool = true) where F <: Function

    check_unique && check_unique_set_name(mesh,set_name,:edge_sets)

    e_set = Set{Int}()
    topo = mesh.topo

    for edge in RootIterator{2}(topo)

        if all(
            f(mesh[node_id]) && is_active(mesh[node_id])
                for node_id in get_edge_node_ids(topo,edge.id)
            )

            push!(e_set,edge.id)
        end
    end

    mesh.edge_sets[set_name] = e_set
end


@kwdef struct ConstraintHandler{D,U,ET}
    mesh::Mesh{D,ET}

    # dirichlet bcs
    d_bcs::Dict{Int,Float64} = Dict{Int,Float64}()

    # neumann bcs
    n_bcs::Dict{Int,Float64} = Dict{Int,Float64}()
end

function ConstraintHandler{U}(mesh::Mesh{D,ET}) where {D,U,ET}
    return ConstraintHandler{D,U,ET}(mesh=mesh)
end



function add_dirichlet_bc!(ch::ConstraintHandler{D,U},
    dh::DofHandler{D,U},
    set_name::String,
    f::F;
    c_dofs::AbstractVector{Int} = SVector{U,Int}(1:U)) where {D,U,F <: Function}


    @assert U >= maximum(c_dofs) "wrong dim of constraint dofs"

    node_set = ch.mesh.node_sets[set_name]

    mesh = ch.mesh

    dummy_input = zero(eltype(mesh.nodes))
    dummy_output = f(dummy_input)

    @assert length(dummy_output) == length(c_dofs) "wrong dim of constraint dofs" 

    for node_id in node_set
        !is_active(mesh[node_id]) && continue
        dof_ids = get_dofs(dh,node_id)[c_dofs]
        
        f_output = f(mesh[node_id]) 

        for (i,dofi) in enumerate(dof_ids)
            ch.d_bcs[dofi] = f_output[i]
        end
    end
    nothing
end


# TODO: add neumann bc for edges and faces


function add_neumann_bc!(ch::ConstraintHandler{D,U},
    set_name::String,
    f::F;
    c_dofs::AbstractVector{Int} = SVector{U,Int}(1:U)) where {D,U,F <: Function}

    @assert U >= maximum(c_dofs) "wrong dim of constraint dofs"

    throw("Not implemented")
end