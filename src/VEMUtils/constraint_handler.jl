function check_unique_set_name(mesh::Mesh,set_name::String,set_type::Symbol)

    set = if set_type == :node_sets
        mesh.node_sets
    elseif set_type == :edge_sets
        mesh.edge_sets
    elseif set_type == :face_sets
        mesh.face_sets
    elseif set_type == :volume_sets
        mesh.volume_sets
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


function get_manifold_set_name(manifold_dim::Val{D}) where D
    if D == 2
        return :edge_sets
    elseif D == 3
        return :face_sets
    elseif D == 4
        return :volume_sets
    end
    error("Invalid manifold dimension")
end

function add_manifold_set!(mesh::Mesh,
    set_name::String, 
    f::F, 
    manifold_dim::Val{D};
    check_unique::Bool = true) where {D,F <: Function}

    manifold_set_name = get_manifold_set_name(manifold_dim)
    check_unique && check_unique_set_name(mesh,set_name,manifold_set_name)

    m_set = Set{Int}() 
    topo = mesh.topo

    for manifold in RootIterator{D}(topo)
        m_node_ids = topo.connectivity[1, D][manifold.id]
        if all(f(mesh[node_id]) && is_active(mesh[node_id]) for node_id in m_node_ids)
            push!(m_set,manifold.id)
        end
    end

    getfield(mesh,manifold_set_name)[set_name] = m_set
end

@inline add_edge_set!(mesh::Mesh,set_name,f::Function;check_unique= true)  = 
                add_manifold_set!(mesh,set_name,f,Val(2);check_unique=check_unique)
@inline add_face_set!(mesh::Mesh,set_name,f::Function;check_unique= true)  = 
                add_manifold_set!(mesh,set_name,f,Val(3);check_unique=check_unique)
@inline add_volume_set!(mesh::Mesh,set_name,f::Function;check_unique=true) = 
                add_manifold_set!(mesh,set_name,f,Val(4);check_unique=check_unique)



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




function function_integration2d(f::Function, 
    points::AbstractVector{<:StaticVector{2,Float64}}; 
    order = 2)

	n_points = length(points)

	alpha = sum(p[1] for p in points)/n_points


    gp,gw = get_gauss_legendre(order+1)
	res = zero(Float64)

	@inbounds for i in eachindex(points)
		next_i = mod1(i+1,n_points)
		ip = points[i]
		np = points[next_i]
		
		for j in 1:order+1
			xj,yj = (np - ip)/2 * gp[j] .+ (np + ip)/2
			wj  = (np[2] - ip[2])/4 * (xj-alpha) * gw[j]
			for m in 1:order+1
				xjm = (xj-alpha)/2 * gp[m] + (xj+alpha)/2
				
				wjm = wj * gw[m]	
				res += f(xjm,yj) * wjm
			end
		end
	end
	return res
end





"""
    add_neumann_bc!(ch::ConstraintHandler{D,U,ET}, dh::DofHandler{D,U}, 
                    fd_col::Dict{Int,<:FaceData{D}}, set_name::String, 
                    f::Function) where {D,U,K,ET<:ElType{K}}

Add Neumann boundary conditions to a constraint handler by integrating a given function 
over faces in a specified face set.

# Arguments
- `ch::ConstraintHandler{D,U,ET}`: The constraint handler to which Neumann BCs will be added
- `dh::DofHandler{D,U}`: The degree of freedom handler for the mesh
- `fd_col::Dict{Int,<:FaceData{D}}`: Dictionary mapping face IDs to face data structures
- `set_name::String`: Name of the face set in the mesh where BCs should be applied
- `f::Function`: Function to integrate over the boundary faces. Should accept 3D coordinates 
                 and return a value compatible with the dot product with base functions

# Details
For each face in the specified face set, this function:
1. Projects the face nodes to 2D coordinates in the face's local plane
2. Retrieves the element base functions for the face
3. Integrates the product of the given function `f` and each base function over the face
4. Accumulates the integrated values into the constraint handler's Neumann BC storage

The integration is performed in the face's local 2D coordinate system, with automatic 
transformation between 3D and 2D coordinates.

# Notes
- Currently only tested for `K == 1`
- The function modifies `ch.n_bcs` in place
- Multiple calls to this function will accumulate values for the same DOFs
"""
function add_neumann_bc!(ch::ConstraintHandler{D,U,ET},
    dh::DofHandler{D,U},
    fd_col::Dict{Int,<:FaceData{D}},
    set_name::String,
    f::F) where {D,U,K,ET<:ElType{K},F <: Function}

    @assert K == 1 "so far only tested for K == 1"
    face_set = ch.mesh.face_sets[set_name] 

    for (face_id,fd) in fd_col 

        face_id ∉ face_set && continue


        face_base      = ElementBaseFunctions(get_base(BaseInfo{2,K,U}()),fd.ΠsL2 |> x -> stretch(x,Val(U)))

        nodes2d = map(fd.face_node_ids) do node_id
            project_to_2d_abs(ch.mesh[node_id],fd.dΩ.plane)
        end


        dofs = get_dofs(dh,fd.face_node_ids)


        for (base_fun,dof) in zip(face_base,dofs)
   
            val = function_integration2d(nodes2d,order = 1) do x,y 
                x2d = SA[x,y]
                x_scaled = (x2d - fd.dΩ.bc) / fd.dΩ.hf 

                x3d = project_to_3d(x2d,fd.dΩ.plane)

                return f(x3d) ⋅ base_fun(x_scaled)
            end
  
            ch.n_bcs[dof] = get(ch.n_bcs, dof, 0.0) + val
        end
    end

    nothing
end