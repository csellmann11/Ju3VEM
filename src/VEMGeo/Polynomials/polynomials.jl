#### Polyomials ##################
include("monomials.jl")
include("../stretched_matrices.jl")


struct Polynomial{T,D,L} <: Number
    coeffs::SVector{L,T}
    base  ::SVector{L,Monomial{Float64,D}}
end

Base.zero(p::Polynomial{T,D,L}) where {T,D,L} = p(zero(SVector{L,T}),p.base)

Base.isequal(p::Polynomial,q::Polynomial) = p.coeffs == q.coeffs && p.base == q.base
Base.isapprox(p::Polynomial,q::Polynomial) = p.coeffs ≈ q.coeffs && p.base == q.base

function Polynomial(coeffs::AbstractVector{T},base::SVector{L,Monomial{Float64,D}}) where {T,D,L}
    return Polynomial{T,D,L}(SVector{L,T}(coeffs),base)
end

@inline (p::Polynomial{T1})(x::V) where {T1<:Real,T2<:Real,V<:AbstractVector{T2}} = sum(c*b(x) for (c,b) in zip(p.coeffs,p.base)) 

@inline function (p::Polynomial{T,D,L})(x::V,bc::StaticVector{D,T},h::T) where {T,D,L,V<:AbstractVector{T}}
    return sum(c*b(x,bc,h) for (c,b) in zip(p.coeffs,p.base)) 
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

@inline get_base_len(D, O, U) = max(div(prod(O+1:O+D), factorial(D)) * U,0)

function get_bi_base(::BaseInfo{D,O,U}) where {D, O, U}
    get_base(BaseInfo{D, O, U}(),Val(true))
end


for D in 1:3, O in 1:6 
    get_base(BaseInfo{D,O,1}())
end
println("done getting base")



@generated function get_order(::Val{L},::Val{D}) where {L,D}
    O = 0 
    while true 
        base_len = get_base_len(D,O,1)
        if base_len == L 
            return O
        else
            O += 1
        end
    end
end

function get_order(p::Polynomial{T,D,L}) where {T,D,L}
    return get_order(Val(L),Val(D))
end

for D in 1:3, O in 1:6
    L = get_base_len(D,O,1)
    get_order(Val(L),Val(D)) 
end
println("done getting order")







# Precompute binomials globally
const MAXN = 5
const MAXD = 20
const BINOMS = [binomial(n, k) for n in 0:MAXN+MAXD, k in 0:MAXD]

@inline function bchoose(n::Int, k::Int)
    (n < 0 || k < 0 || k > n) && return 0
    return BINOMS[n+1, k+1]
end
# prefix sums: number of monomials of degree ≤ d
const BINOMSUMS = [sum(binomial(n + k - 1, k) for k in 0:d) for n in 1:MAXN, d in 0:MAXD]

@inline function get_exp_to_idx_dict(_exp::SVector{N,Int},::Val{false} = Val(false)) where {N}
    # @assert N <= MAXN "increase MAXN"
    exp = reverse(_exp)             # keep same reversal convention as you used
    d = sum(exp)
    # @assert d <= MAXD "increase MAXD"

    # count all monomials of strictly smaller total degree
    idx = (d == 0) ? 0 : BINOMSUMS[N, d]  # BINOMSUMS[N,d] stores sum_{k=0}^{d-1} C(N+k-1,k)
                                          # note: we stored BINOMSUMS with second index = d+1 -> adjusted below
    # The above line uses BINOMSUMS[N, d] because we filled BINOMSUMS[N, d+1] = sum_{k=0}^d ...
    # So BINOMSUMS[N, d] == sum_{k=0}^{d-1} ...
    used = 0
    @inbounds for i in 1:N
        a = exp[i]
        if a > 0
            S = d - used
            # closed form: sum_{t=0}^{a-1} C(N-i + S - t - 1, S - t)
            # = C(N - i + S, S) - C(N - i + S - a, S - a)
            n1 = N - i + S        # top for term1
            k1 = S
            n2 = N - i + S - a    # top for term2
            k2 = S - a
            # Map to BINOMS: C(n,k) == BINOMS[n+1, k+1]
            term1 = BINOMS[n1 + 1, k1 + 1]
            term2 = BINOMS[n2 + 1, k2 + 1]
            idx += term1 - term2
        end
        used += a
    end

    return idx + 1
end





function Base.:*(p::Polynomial{T1,D,L1},q::Polynomial{T1,D,L2}) where {T1,D,L1,L2}
    O = get_order(p) + get_order(q)
    new_base = get_base(BaseInfo{D,O,1}()).base
    # _coeffs = zero(MVector{length(new_base),T1})
    _coeffs = FixedSizeVector{T1}(undef,length(new_base))
    _coeffs .= zero(T1)
    for (i,mi) in enumerate(p.base)
        ci = p.coeffs[i]
        # ci == zero(ci) & continue
        for (j,mj) in enumerate(q.base)
            m = mi * mj
            idx = get_exp_to_idx_dict(m.exp)
            _coeffs[idx] += ci * q.coeffs[j]
        end
    end
    coeffs = SVector{length(new_base)}(_coeffs)
    return Polynomial(coeffs,new_base)
end

function poly_pow(p::Polynomial{T,D,L},::Val{n}) where {T,D,L,n} 
    @assert n >= 0 "Exponent must be non-negative; got $(n)"
    n == 1 && return p 
    # n == 0 && return Polynomial(SVector{L}(i == 1 ? one(T) : zero(T) for i in 1:L),p.base)
    n == 0 && return Polynomial(SA[one(T)],get_base(BaseInfo{D,0,1}()).base)

    if isodd(n)
        return p * poly_pow(p,Val(n-1))
    else
        p = poly_pow(p,Val(n÷2))
        return p * p
    end
end





function extract_constructor(::Type{T}) where T<:AbstractArray
    return T.name.wrapper
end
const MonoPoly{T,D} = Union{Polynomial{T,D},Monomial{T,D}}

@inline Base.Vector(itr::Base.Generator) = collect(itr)

function (v::AbstractVector{<:MonoPoly})(x)
    U = length(typeof(v)) # does only work for statically sized vectors 
    return extract_constructor(typeof(v))(SVector{U}(v[i](x) for i in 1:U))
end

function (v::Union{Vector{M},FixedSizeVector{M}} where M<:MonoPoly)(x) 
    # return extract_constructor(typeof(v))(v[i](x) for i in 1:length(v))
    # collect_as(extract_constructor(typeof(v)),v[i](x) for i in eachindex(v))
    out = FixedSizeVector{Float64}(undef,length(v))
    for i in eachindex(v)
        out[i] = v[i](x)
    end
    return out
end


function (v::AbstractVector{<:MonoPoly})(x,bc,h)
    U = length(typeof(v)) # does only work for statically sized vectors 
    return extract_constructor(typeof(v))(SVector{U}(v[i](x,bc,h) for i in 1:U))
end

function (v::AbstractArray{<:MonoPoly})(x)
    return extract_constructor(typeof(v))(SArray{Tuple{size(v)...}}(v[i](x) for i in eachindex(v)))
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






function poly_grad(p::Polynomial{T,D,L},h::T) where {T,D,L}
    O = get_order(p)
    grad_base  = get_base(BaseInfo{D,O-1,1}())
    grad_coeffs = SVector{2}(zero(MVector{length(grad_base),T}) for _ in 1:2)

    @inbounds for i in 2:L 
        nabm = ∇(p.base[i],h) 
        for d in 1:D 
            exp = getexp(nabm[d])
            val = getval(nabm[d]) 
            val == zero(T) && continue
            idx = get_exp_to_idx_dict(exp)
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
