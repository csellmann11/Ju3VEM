#### Polyomials ##################
include("monomials.jl")
include("../stretched_matrices.jl")


@inline get_base_len(D, O, U) = max(div(prod(O+1:O+D), factorial(D)) * U, 0)

@generated function get_order(::Val{L}, ::Val{D}) where {L,D}
    O = 0
    while true
        base_len = get_base_len(D, O, 1)
        if base_len == L
            return O
        else
            O += 1
        end
    end
end




# this is needed to avoid world age issues
for D in 1:3, O in 1:6
    L = get_base_len(D, O, 1)
    get_order(Val(L), Val(D))
end


struct Polynomial{T,D,L} <: Number
    coeffs::SVector{L,T}
    base::SVector{L,Monomial{Float64,D}}
end

Base.zero(p::Polynomial{T,D,L}) where {T,D,L} = Polynomial(zero(SVector{L,T}), p.base)

Base.isequal(p::Polynomial, q::Polynomial) = p.coeffs == q.coeffs && p.base == q.base
Base.isapprox(p::Polynomial, q::Polynomial) = p.coeffs ≈ q.coeffs && p.base == q.base

function Polynomial(coeffs::AbstractVector{T}, base::SVector{L,Monomial{Float64,D}}) where {T,D,L}
    return Polynomial{T,D,L}(SVector{L,T}(coeffs), base)
end


function get_order(::Type{Polynomial{T,D,L}}) where {T,D,L}
    return get_order(Val(L), Val(D))
end

function get_order(p::Polynomial{T,D,L}) where {T,D,L}
    return get_order(typeof(p))
end



@generated function poly_eval(p::Polynomial{T,D,L},
    x::V,
    bc::StaticVector{D,T},
    h::T) where {T,D,L,V<:AbstractVector{T}}
    O = get_order(p) 
    base = get_base(BaseInfo{D,O,1}())
    
    expr = Expr(:call, :+, 0)
    for i in 1:L
        m = base[i]
        inner_expr = Expr(:call, :*, :(p.coeffs[$i]))
        for d in 1:D
            m.exp[d] == 0 && continue
            push!(inner_expr.args, :(x_scaled[$d]^$(m.exp[d])))  # Use x_scaled!
        end
        push!(expr.args, inner_expr)
    end
    
    return quote
        x_scaled = (x .- bc) / h
        $expr
    end
end


# @inline (p::Polynomial{T1})(x::V) where {T1<:Real,T2<:Real,V<:AbstractVector{T2}} = sum(c * b(x) for (c, b) in zip(p.coeffs, p.base))
@inline (p::Polynomial{T1})(x::V) where {T1<:Real,T2<:Real,V<:AbstractVector{T2}} = poly_eval(p, x, zero(x), 1.0)

@inline function (p::Polynomial{T,D,L})(x::V, bc::StaticVector{D,T}, h::T) where {T,D,L,V<:AbstractVector{T}}
    # return sum(c * b(x, bc, h) for (c, b) in zip(p.coeffs, p.base))
    return poly_eval(p, x, bc, h)
end











struct PolynomialBase{D,O,U,L}
    base::SVector{L,Monomial{Float64,D}}
end
function PolynomialBase{O,U}(base::SVector{L,Monomial{Float64,D}}) where {D,O,U,L}
    PolynomialBase{D,O,U,L}(base)
end
function Base.getindex(pbase::PolynomialBase{D,O,U}, idx::Int) where {D,O,U}
    U == 1 && return pbase.base[idx]

    base_idx, v_idx = divrem(idx - 1, U) .+ (1, 1)
    return SVector{U}(i != v_idx ? zero(Monomial{Float64,D}) : pbase.base[base_idx] for i in 1:U)
end

Base.length(pb::PolynomialBase{D,O,U}) where {D,O,U} = length(pb.base) * U
Base.eltype(pb::PolynomialBase{D,O,1}) where {D,O} = typeof(pb.base[1])
Base.eltype(pb::PolynomialBase{D,O,U}) where {D,O,U} = SVector{U,typeof(pb.base[1])}

function Base.iterate(pb::PolynomialBase{D,O,U}, state::Int=1) where {D,O,U}
    if state > length(pb)
        return nothing
    else
        return (pb[state], state + 1)
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
@generated function get_base(::BaseInfo{D,O,U}, ::Val{BB}=Val(false)) where {D,O,U,BB}
    @assert isa(BB, Bool) "SB must be bool or true"

    base_len = if BB
        (O + 1)^D
    else
        div(prod(O+1:O+D), factorial(D))
    end


    base = Vector{Monomial{Float64,D}}(undef, base_len)
    count = 1
    max_order = O + BB
    for order in 0:max_order, n in CartesianIndices(ntuple(_ -> order + 1, Val(D)))
        exp = SVector(n.I .- 1)
        (sum(exp) != order || max(exp...) > O) && continue
        base[count] = Monomial(1.0, exp)
        count += 1
    end
    base = PolynomialBase{O,U}(SVector{base_len}(base))

    return :($base)
end



function get_bi_base(::BaseInfo{D,O,U}) where {D,O,U}
    get_base(BaseInfo{D,O,U}(), Val(true))
end

# this is needed to avoid world age issues
for D in 1:3, O in 1:6
    get_base(BaseInfo{D,O,1}())
end










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

@inline function get_exp_to_idx_dict(_exp::SVector{N,Int}, ::Val{false}=Val(false)) where {N}
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
            term1 = BINOMS[n1+1, k1+1]
            term2 = BINOMS[n2+1, k2+1]
            idx += term1 - term2
        end
        used += a
    end

    return idx + 1
end






function Base.:*(p::Polynomial{T1,D,L1}, q::Polynomial{T1,D,L2}) where {T1,D,L1,L2}
    O = get_order(p) + get_order(q)
    new_base = get_base(BaseInfo{D,O,1}()).base
    # _coeffs = zero(MVector{length(new_base),T1})
    _coeffs = FixedSizeVector{T1}(undef, length(new_base))
    _coeffs .= zero(T1)
    for (i, mi) in enumerate(p.base)
        ci = p.coeffs[i]
        # ci == zero(ci) & continue
        for (j, mj) in enumerate(q.base)
            m = mi * mj
            idx = get_exp_to_idx_dict(m.exp)
            _coeffs[idx] += ci * q.coeffs[j]
        end
    end
    coeffs = SVector{length(new_base)}(_coeffs)
    return Polynomial(coeffs, new_base)
end

function poly_pow(p::Polynomial{T,D,L}, ::Val{n}) where {T,D,L,n}
    @assert n >= 0 "Exponent must be non-negative; got $(n)"
    n == 1 && return p
    # n == 0 && return Polynomial(SVector{L}(i == 1 ? one(T) : zero(T) for i in 1:L),p.base)
    n == 0 && return Polynomial(SA[one(T)], get_base(BaseInfo{D,0,1}()).base)

    if isodd(n)
        return p * poly_pow(p, Val(n - 1))
    else
        p = poly_pow(p, Val(n ÷ 2))
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
    out = FixedSizeVector{Float64}(undef, length(v))
    for i in eachindex(v)
        out[i] = v[i](x)
    end
    return out
end


function (v::AbstractVector{<:MonoPoly})(x, bc, h)
    U = length(typeof(v)) # does only work for statically sized vectors 
    return extract_constructor(typeof(v))(SVector{U}(v[i](x, bc, h) for i in 1:U))
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
    sol::AbstractVector, Π_star::StretchedMatrix) where {D,O,U,L}


    v = SVector{U}(begin
        c = matmul(Π_star.data, @view(sol[i:U:end]))
        Polynomial(c, b.base)
    end for i in 1:U)

    if U == 1
        return only(v)
    else
        return v
    end
end

function idx_rem_div(v, U)
    q = (v + U - 1) ÷ U
    r = mod1(v, U)   
    return q, r 
end

function get_static_col(M::AbstractMatrix,
    col_idx::Int,
    ::Val{L}) where {L}

    @assert size(M, 1) == L "M must have L rows, got $(size(M,1))"
    @assert 0 < col_idx <= size(M, 2) "col_idx must be between 1 and $(size(M,2)), got $col_idx"
    out = zero(MVector{L,Float64})
    for i in 1:L
        out[i] = M[i, col_idx]
    end
    return SVector(out)
end

function get_static_row(M::AbstractMatrix,
    row_idx::Int, ::Val{L}) where {L}
    @assert size(M, 2) == L "M must have L columns, got $(size(M,2))"
    @assert 0 < row_idx <= size(M, 1) "row_idx must be between 1 and $(size(M,1)), got $row_idx"
    out = zero(MVector{L,Float64})
    for i in 1:L
        out[i] = M[row_idx, i]
    end
    return SVector(out)
end

function unit_sol_proj(b::PolynomialBase{D,O,U,L},
    one_idx::Int, Π_star::StretchedMatrix) where {D,O,U,L}


    proj_col, v_idx = idx_rem_div(one_idx, U)
    coeffs          = get_static_col(Π_star.data, proj_col, Val(L))

    p               = Polynomial(coeffs, b.base)

    if U == 1
        return p
    else
        return SVector{U}(i == v_idx ? p : zero(p) for i in 1:U)
    end
end


function grad_sol_proj(b::PolynomialBase{D,O,U,L},
    sol::AbstractVector, Π_star::StretchedMatrix) where {D,O,U,L}

    @assert U == D "U must equal D, got D = $D and U = $U"
    len = length(sol)
    v = Mat{U,D}((i, j) -> begin
        c = SVector{2 * length(b.base)}(matmul(Π_star.data, view(sol, i:U:len)))[j:D:end]
        Polynomial(c, b.base)
    end)
    v
end

function get_ϕi(b::PolynomialBase{D,O,U,L}, Π_star::StretchedMatrix{T}, idx::Int) where {D,O,U,L,T}

    base_idx, v_idx = divrem(idx - 1, U) .+ (1, 1)

    coeffs = SVector{L,T}(Π_star.data[j, base_idx] for j in 1:L)

    U == 1 && return Polynomial(coeffs, b.base)

    zero_p = Polynomial(zero(SVector{L,Float64}), b.base)
    return Vec{U}(i -> i == v_idx ? Polynomial(coeffs, b.base) : zero_p)
end






# function poly_grad(p::Polynomial{T,D,L}, h::T) where {T,D,L}
#     O = get_order(p)
#     grad_base = get_base(BaseInfo{D,O - 1,1}())
#     vL = Val(length(grad_base))

#     @no_escape begin
#         grad_coeffs = @alloc(T,D,length(grad_base))
#         grad_coeffs .= zero(T)
     
#         for i in 2:L 
#             m = p.coeffs[i] * p.base[i]
#             ∇m = ∇(m, h)
#             for d in 1:D
#                 exp = getexp(∇m[d])
#                 val = getval(∇m[d])
#                 val == zero(T) && continue
#                 idx = get_exp_to_idx_dict(exp)
#                 grad_coeffs[d,idx] += val
#             end
#         end
#         SVector{D}(Polynomial(get_static_row(grad_coeffs,i,vL), grad_base.base) for i in 1:D)
#     end
# end

@generated function poly_grad(p::Polynomial{T,D,L}, h::T) where {T,D,L}
    O = get_order(p)
    base = get_base(BaseInfo{D,O,1}())
    grad_base = get_base(BaseInfo{D,O-1,1}()).base
    grad_L = length(grad_base)
    
    # Precompute contribution map at compile time
    # contributions[d][grad_idx] = [(poly_idx, factor), ...]
    contributions = [Vector{Tuple{Int,typeof(base[1].val)}}[] for _ in 1:D]
    for _ in 1:grad_L
        for d in 1:D
            push!(contributions[d], Tuple{Int,typeof(base[1].val)}[])
        end
    end
    
    for i in 2:L  # Skip constant term (index 1)
        m = base[i]
        for d in 1:D
            m.exp[d] == 0 && continue
            
            # Compute new exponent after differentiation
            # new_exp = copy(m.exp)
            # new_exp[d] -= 1
            new_exp = setindex(m.exp, m.exp[d] - 1, d)
            coeff_factor = m.val * m.exp[d]
            
            # Find index in gradient base
            grad_idx = findfirst(gm -> all(gm.exp .== new_exp), grad_base)
            grad_idx === nothing && continue
            
            push!(contributions[d][grad_idx], (i, coeff_factor))
        end
    end
    
    # Generate code for each gradient component
    grad_poly_exprs = []
    for d in 1:D
        coeff_exprs = []
        for grad_idx in 1:grad_L
            if isempty(contributions[d][grad_idx])
                push!(coeff_exprs, :(zero(T)))
            else
                terms = [:(p.coeffs[$i] * $(factor)) for (i, factor) in contributions[d][grad_idx]]
                push!(coeff_exprs, Expr(:call, :+, terms...))
            end
        end
        push!(grad_poly_exprs, :(Polynomial(SVector{$grad_L,T}($(coeff_exprs...)) / h, grad_base)))
    end
    
    return quote
        grad_base = get_base(BaseInfo{$D,$(O-1),1}()).base
        SVector{$D}($(grad_poly_exprs...))
    end
end

function poly_grad(p::V, h::T) where {T,V<:AbstractArray{<:Polynomial{T,D,L}}} where {D,L}
    v_size = size(typeof(p)) 
    _grad = SArray{Tuple{v_size...}}(poly_grad(p[i], h) for i in eachindex(p))

    new_v_size = (v_size..., D)
    return SArray{Tuple{new_v_size...}}(begin
        (v1..., v2) = idx.I
        _grad[CartesianIndex(v1)][v2]
    end for idx in CartesianIndices(new_v_size))
end



const ∇p = poly_grad








function gradient_x(p::M, h::T, x::V) where {M<:Union{AbstractVector{<:MonoPoly},MonoPoly},T,V<:AbstractVector{T}}
    ∇p(p, h)(x)
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
