

mutable struct CellValues{D,U,ET,L1,L2,V<:AbstractVector}

    const mesh::Mesh{D,ET}
    const vnm ::NodeID2LocalID

    const facedata_col::Dict{Int,FaceData{D,L1,V}}
    volume_data ::VolumeIntegralData{L2}

    const dh::DofHandler{D,U}
end



function _create_facedata_col(mesh::Mesh{D,StandardEl{K}}) where {D,K}
    K_MAX = max(2*K-1, 2)
    topo = mesh.topo


    FI_type = Core.Compiler.return_type(precompute_face_monomials,
            Tuple{Int,Topology{D},Val{K_MAX}})

    FD_Type = Core.Compiler.return_type(h1_projectors!,
            Tuple{Int,Mesh{3,StandardEl{K}},FI_type})
    facedata_col = Dict{Int,FD_Type}()

    for face in RootIterator{3}(topo)
        facedata_col[face.id] = h1_projectors!(
            face.id, 
            mesh, 
            precompute_face_monomials(face.id, topo, Val(K_MAX))
            )
    end
    return facedata_col
end


function CellValues{U}(
    mesh::Mesh{D,ET}) where {D,U,K,ET<:ElType{K}}

    dh = DofHandler{U}(mesh) 

    vnm = NodeID2LocalID(;sizehint = 30) 

    # base is just a dummy 3d base, even if problem is 2d
    base = get_base(BaseInfo{3,K,1}()).base
    volume_data = VolumeIntegralData{length(base)}(
        0.0,
        zero(SVector{3,Float64}),
        zero(SVector{length(base),Float64}))

    facedata_col = _create_facedata_col(mesh)
    V = typeof(first(values(facedata_col)).face_node_ids)

    return CellValues(mesh,vnm,facedata_col,volume_data,dh)
end

function get_element_bc(cv::CellValues)
    return cv.volume_data.vol_bc
end

function get_element_diameter(cv::CellValues)
    return cv.volume_data.hvol
end

function get_element_area(cv::CellValues)
    return cv.volume_data.integrals[1]
end


function reinit!(el_id::Int,cv::CellValues{3})
    
    mesh = cv.mesh; topo = mesh.topo
    K = get_order(cv.mesh)

    vol_data = precompute_volume_monomials(
        el_id, topo, cv.facedata_col, Val(K))

    create_node_mapping!(cv.vnm,el_id,mesh,cv.facedata_col)

    cv.volume_data = vol_data

    return cv
end

@inline get_n_cell_dofs(cv::CellValues{D,U}) where {D,U} = U*get_n_nodes(cv.vnm)

function get_cell_dofs(cv::CellValues{D,U}) where {D,U}
    dofs = FixedSizeVector{Int}(undef, U*get_n_nodes(cv.vnm))
    for (node_id,i) in cv.vnm.map
        node_dofs = get_dofs(cv.dh,node_id)
        for j in 1:U
            idx = (i-1)*U + j
            dofs[idx] = node_dofs[j]
        end
    end
    return dofs
end




"""
    ElementBaseFunctions{D,O,U,L,M,V} <: AbstractVector{V}

A vector-like container representing the projected element base functions for a VEM element.

# Type Parameters
- `D`: Spatial dimension of the element
- `O`: Order of the polynomial basis
- `U`: Number of unknowns per DOF (1 for scalar problems, >1 for vector-valued problems)
- `L`: Length of the polynomial basis (derived from D and O)
- `M`: Type of the projection matrix (must be a subtype of `StretchedMatrix`)
- `V`: Element type, either `Polynomial{Float64,D,L}` for scalar problems or 
       `SVector{U,Polynomial{Float64,D,L}}` for vector-valued problems

# Fields
- `base::PolynomialBase{D,O,U,L}`: The polynomial basis for the element
- `Π_star::M`: The projection matrix (Π*) that maps DOF values to polynomial coefficients

# Description
This struct provides an abstract vector interface to access individual element base functions.
Each base function is obtained by projecting a unit DOF Vector using the `Π_star` matrix.
For scalar problems (U=1), each element is a single polynomial. For vector-valued problems 
(U>1), each element is an SVector of polynomials.

# Interface
- `length(ebf)`: Returns the number of base functions (DOFs) for the element
- `ebf[i]`: Returns the i-th projected base function as a Polynomial or SVector of Polynomials
- Supports iteration through all base functions

# Example Usage
```julia
# Access the i-th base function
ϕ_i = element_base_functions[i]

# Evaluate at a point
value = ϕ_i(x)  # where x is a coordinate in the element's local coordinate system
```

# Constructor
"""
struct ElementBaseFunctions{D,
    O,U,L,
    M<:StretchedMatrix} <: AbstractVector{SVector{U,Polynomial{Float64,D,L}}}

    base::PolynomialBase{D,O,U,L}
    Π_star::M


    function ElementBaseFunctions(base::PolynomialBase{D,O,U,L}, Π_star::M) where {D,O,U,L,M}
        new{D,O,U,L,M}(base, Π_star)
    end
end

Base.length(ebf::ElementBaseFunctions) = size(ebf.Π_star, 2)
Base.size(ebf::ElementBaseFunctions) = (size(ebf.Π_star, 2),)
Base.eltype(::Type{ElementBaseFunctions{D,O,U,L,M}}) where {D,O,U,L,M} = SVector{U,Polynomial{Float64,D,L}}
Base.getindex(ebf::ElementBaseFunctions, idx::Int) = unit_sol_proj(ebf.base,idx,ebf.Π_star)


# function get_base_fun(ebf::ElementBaseFunctions, 
#     idx::Int)
#     return unit_sol_proj(ebf.base,idx,ebf.Π_star)
# end

# # iteration interface 
# function Base.iterate(ebf::ElementBaseFunctions, state::Int=1)
#     N = size(ebf.Π_star, 2) 
#     state > N && return nothing
#     return (get_base_fun(ebf, state), state + 1)
# end

# Base.length(ebf::ElementBaseFunctions) = size(ebf.Π_star, 2)

# Base.eltype(::Type{ElementBaseFunctions{D,O,1,L,M}}) where {D,O,L,M} =
#     Polynomial{Float64,D,L}  

# Base.eltype(::Type{ElementBaseFunctions{D,O,U,L,M}}) where {D,O,U,L,M} =
#     SVector{U,Polynomial{Float64,D,L}}