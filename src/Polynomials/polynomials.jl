#### Polyomials ##################
include("monomials.jl")
include("../stretched_matrices.jl")


struct Polynomial{T,D,L} <: Number
    coeffs::SVector{L,T}
    base  ::SVector{L,Monomial{Float64,D}}
end

Base.zero(p::Polynomial{T,D,L}) where {T,D,L} = p(zero(SVector{L,T}),p.base)

function Polynomial(coeffs::AbstractVector{T},base::SVector{L,Monomial{T,D}}) where {T,D,L}
    return Polynomial{T,D,L}(SVector{L,T}(coeffs),base)
end

@inline (p::Polynomial{T1})(x::V) where {T1<:Real,T2<:Real,V<:AbstractVector{T2}} = sum(c*b(x) for (c,b) in zip(p.coeffs,p.base)) 



const EXP_TO_IDEX_DICT_D1 = Dict{SVector{1,Int},Int}()
const EXP_TO_IDEX_DICT_D2 = Dict{SVector{2,Int},Int}()
const EXP_TO_IDEX_DICT_D3 = Dict{SVector{3,Int},Int}()

function get_exp_to_idx_dict(exp::SVector{N,Int},O::Int) where N
    exp_dict = if N == 1
        EXP_TO_IDEX_DICT_D1
    elseif N == 2
        EXP_TO_IDEX_DICT_D2
    elseif N == 3
        EXP_TO_IDEX_DICT_D3
    else
        throw(ErrorException("Only supported for N == 1, 2 or 3"))
    end

    get!(exp_dict,exp) do 
        base = get_base(BaseInfo{N,O,1}(),Val(false))
        for (i,m) in enumerate(base.base)
            if m.exp == exp 
                return i
            end
        end
        throw(ErrorException("Exponent $exp not found in base"))        
    end
end



function extract_constructor(::Type{T}) where T<:AbstractArray
    return T.name.wrapper
end
const MonoPoly{T,D} = Union{Polynomial{T,D},Monomial{T,D}}

function (v::AbstractVector{<:MonoPoly})(x)
    U = length(typeof(v)) # does only work for statically sized vectors 
    return extract_constructor(typeof(v))(SVector{U}(v[i](x) for i in 1:U))
end



function (v::AbstractVector{<:MonoPoly})(x,bc,h)
    U = length(typeof(v)) # does only work for statically sized vectors 
    return extract_constructor(typeof(v))(SVector{U}(v[i](x,bc,h) for i in 1:U))
end

function (v::AbstractArray{<:MonoPoly})(x)
    return extract_constructor(typeof(v))(SArray{Tuple{size(v)...}}(v[i](x) for i in eachindex(v)))
end




struct PolynomialBase{D,O,U,L} 
    base::SVector{L,Monomial{Float64,D}}
end
function PolynomialBase{O,U}(base::SVector{L,Monomial{Float64,D}}) where {D,O,U,L}
    PolynomialBase{D,O,U,L}(base)
end
function Base.getindex(pbase::PolynomialBase{D,O,U},idx::Int) where {D,O,U}
    U == 1 && return pbase.base[idx]

    base_idx, v_idx = divrem(idx-1,U) .+ (1,1)
    return SVector{U}(i != v_idx ? zero(Monomial{Float64,D}) : pbase.base[base_idx] for i in 1:U)
end

Base.length(pb::PolynomialBase{D,O,U}) where {D,O,U} = length(pb.base)*U
Base.eltype(pb::PolynomialBase{D,O,1}) where {D,O}   = typeof(pb.base[1])
Base.eltype(pb::PolynomialBase{D,O,U}) where {D,O,U} = SVector{U,typeof(pb.base[1])}

function Base.iterate(pb::PolynomialBase{D,O,U}, state::Int = 1) where {D,O,U}
    if state > length(pb)
        return nothing
    else
        return (pb[state],state+1)
    end
end


function (pbase::PolynomialBase)(x::V) where {V<:AbstractVector} 
    return sum(pbase[i](x) for i in eachindex(pbase))
end


# Function to generate the polynomial base
"""
    get_base(::BaseInfo{D, O, U},::Val{BB} = Val(false) where {D, O, U} -> PolynomialBase{D,O,U,L}

Generate the polynomial base for given dimensions `D`, order `O`, and unknowns `U`.

# Arguments
- `bi::BaseInfo{D, O, U}`: Base information.
- `::Val{BB} = Val(false)`: wheather to build a bilinear base or not

# Returns
- `base::PolynomialBase{D,O,U,L}` - Polynomial base.
"""
@generated function get_base(::BaseInfo{D, O, U},::Val{BB} = Val(false)) where {D, O, U, BB}
    @assert isa(BB,Bool) "SB must be bool or true"

    base_len = if BB 
        (O+1)^D
    else 
        div(prod(O+1:O+D), factorial(D)) 
    end

    
    base = Vector{Monomial{Float64, D}}(undef, base_len)
    count = 1
    max_order = O + BB
    for order in 0:max_order, n in CartesianIndices(ntuple(_ -> order+1,Val(D)))
        exp = SVector(n.I .- 1)
        (sum(exp) != order || max(exp...) > O) && continue
        base[count] = Monomial(1.0, exp)
        count += 1
    end
    base = PolynomialBase{O,U}(SVector{base_len}(base))

    return :($base)
end


function get_bi_base(::BaseInfo{D,O,U}) where {D, O, U}
    get_base(BaseInfo{D, O, U}(),Val(true))
end




"""
    sol_proj(base::StaticVector{L},sol::AbstractVector,Π_star::AbstractMatrix) where L

Function computes the solution into the polynomial space by \\
`base ⊙ (Π_star * sol)`

# Arguments
- `base::StaticVector{L}`: Polynomial base.
- `sol::AbstractVector`: Solution vector.
- `Π_star::AbstractMatrix`: Projection matrix.

# Returns
- `sol_proj::SVector{P}`: Projected solution where typeof(P) = typeof(sum(base)).
"""
function sol_proj(b::PolynomialBase{D,O,U,L},
    sol::AbstractVector,Π_star::StretchedMatrix) where {D,O,U,L}


    v = SVector{U}(begin 
        c = matmul(Π_star.data,@view(sol[i:U:end]))
        Polynomial(c,b.base)
    end for i in 1:U)

    if U == 1 return only(v)  
    else return v end
end 

function grad_sol_proj(b::PolynomialBase{D,O,U,L},
    sol::AbstractVector,Π_star::StretchedMatrix) where {D,O,U,L}

    @assert U == D "U must equal D, got D = $D and U = $U"
    len = length(sol)
    v = Mat{U,D}((i,j)-> begin 
        c = SVector{2*length(b.base)}(matmul(Π_star.data,view(sol,i:U:len)))[j:D:end]
        Polynomial(c,b.base)
    end)
    v
end 

function get_ϕi(b::PolynomialBase{D,O,U,L},Π_star::StretchedMatrix{T},idx::Int) where {D,O,U,L,T}

    base_idx, v_idx = divrem(idx-1,U) .+ (1,1)

    coeffs = SVector{L,T}(Π_star.data[j,base_idx] for j in 1:L)

    U == 1 && return Polynomial(coeffs,b.base)

    zero_p = Polynomial(zero(SVector{L,Float64}),b.base)
    return Vec{U}(i -> i == v_idx ? Polynomial(coeffs,b.base) : zero_p)
end 



@inline get_base_len(D, O, U) = div(prod(O+1:O+D), factorial(D)) * U
@generated function get_order(::Polynomial{T,D,L}) where {T,D,L}
    O = 0 
    while true 
        base = get_base(BaseInfo{D,O,1}(),Val(false))
        if length(base) == L 
            return O
        else
            O += 1
        end
    end
end


function poly_grad(p::Polynomial{T,D,L},h::T) where {T,D,L}
    O = get_order(p)
    grad_base = get_base(BaseInfo{D,O-1,1}())
    grad_coeffs = SVector{2}(zero(MVector{length(grad_base),T}) for _ in 1:2)

    for i in 2:L 
        nabm = ∇(p.base[i],h) 
        for d in 1:D 
            exp = getexp(nabm[d])
            val = getval(nabm[d]) 
            val == zero(T) && continue
            idx = get_exp_to_idx_dict(exp,O)
            grad_coeffs[d][idx] += val * p.coeffs[i]
        end
    end
    return SVector{D}(Polynomial(grad_coeffs[i],grad_base.base) for i in 1:D)
end

function poly_grad(p::V,h::T) where {T,V<:AbstractArray{<:Polynomial{T,D,L}}} where {D,L}
    v_size = size(typeof(p))
    _grad = SArray{Tuple{v_size...}}(poly_grad(p[i],h) for i in eachindex(p))

    new_v_size = (v_size...,D)
    return SArray{Tuple{new_v_size...}}(begin 
        (v1...,v2) = idx.I
        _grad[CartesianIndex(v1)][v2]
    end for idx in CartesianIndices(new_v_size))
end



const ∇p = poly_grad








function gradient_x(p::M,h::T, x::V) where {M<:Union{AbstractVector{<:MonoPoly},MonoPoly},T,V<:AbstractVector{T}}
    grad = ∇p(p,h)(x)
    return grad/h
end

const ∇x = gradient_x



########## Display Functionality ##############

# String representation for polynomials, leveraging monomial formatting
function _polynomial_string(p::Polynomial{T,D,L}) where {T,D,L}
    result = ""
    first_term = true

    for i in 1:L
        c = p.coeffs[i]
        c == 0 && continue

        exp = p.base[i].exp
        term_str = _monomial_string(Monomial(abs(c), exp))

        if first_term
            if c < 0
                result *= "-" * term_str
            else
                result *= term_str
            end
            first_term = false
        else
            if c < 0
                result *= " - " * term_str
            else
                result *= " + " * term_str
            end
        end
    end

    return result == "" ? "0" : result
end

function Base.show(io::IO, p::Polynomial{T,D,L}) where {T,D,L}
    print(io, _polynomial_string(p))
end

function Base.show(io::IO, ::MIME"text/plain", p::Polynomial{T,D,L}) where {T,D,L}
    println(io, _polynomial_string(p))
end
