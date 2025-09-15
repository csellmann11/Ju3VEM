using StaticArrays

# Abstract type to represent lazy evaluation of mathematical expressions
struct BaseInfo{D, O, U} end

"""
    struct Monomial{T,D} <: LazyType{T,D,1}
        val::T
        exp::SVector{D,Int}
    end

Structure representing a monomial (a single term in a polynomial).
"""
# Structure representing a monomial (a single term in a polynomial)
struct Monomial{T,D} <: Number
    val::T
    exp::SVector{D,Int}
end 
 

getexp(m::Monomial) = m.exp
getval(m::Monomial) = m.val

# getexp(v::StaticVec{U,Monomial{T,D}}) where {U,T,D} = Vec{U}(i->getexp(v[i]))
# getval(v::StaticVec{U,Monomial{T,D}}) where {U,T,D} = Vec{U}(i->getval(v[i]))

Base.zero(::Type{Monomial{T,D}}) where {T,D} = 
        Monomial(zero(T), SVector{D}(zero(Int) for _ in 1:D))
Base.zero(::Monomial{T,D}) where {T,D} = zero(Monomial{T,D})

Base.one(::Type{Monomial{T,D}}) where {T,D} = Monomial(one(T),SVector{D}(zero(Int) for _ in 1:D))
Base.one(::Monomial{T,D}) where {T,D} = one(Monomial{T,D})

Base.conj(m::Monomial{T,D}) where {T,D} = Monomial(conj(m.val),m.exp)


########## Basic Operations ##############
"""
    (m::Monomial{T,D})(x::StaticVector{D,T}) where {T<:Real,D}

Evaluate a monomial at a given point `x`.
"""
function (m::Monomial{T1,D})(x::V) where {T1<:Real,T2<:Real,D,V<:AbstractVector{T2}}
    @assert length(x) == D "length of x must be equal to D, got $(length(x))"
    prod(xi^expi for (xi,expi) in zip(x,m.exp); init = m.val)
end

"""
    (m::Monomial{T1,D})(x::V1,bc::V2,h::T) where {T<:Real,T1<:Real,T2<:Real,T3<:Real,D,
    V1<:AbstractVector{T2},V2<:AbstractVector{T3}}

Evaluate a monomial at the scaled location `(x-bc)/h`

# Arguments
- `m::Monomial{T1,D}`: Monomial to evaluate.
- `x::V1`: Point to evaluate at.
- `bc::V2`: Offset
- `h::T`: Scale factor.
"""
function (m::Monomial{T1,D})(x::V1,bc::V2,h::T) where {T<:Real,T1<:Real,T2<:Real,T3<:Real,D,
    V1<:AbstractVector{T2},V2<:AbstractVector{T3}}

    @assert length(x) == D "length of x must be equal to D, got $(length(x))"
    @assert length(bc) == D "length of bc must be equal to D, got $(length(bc))"

    prod(((xi-bci)/h)^expi for (xi,expi,bci) in zip(x,m.exp,bc); init = m.val)
end

import Base: heads,tail

function derivative(exp::SVector{D,Int},val::T,h::S,idxs::Vararg{Int,N}) where {T,D,S,N}
    idx = heads(idxs)[1]
    @boundscheck @assert idx <= D "index out of bounds"
    exp_val = exp[idx]
    new_exp_val = exp_val-1
    new_exp = if new_exp_val < 0
        zero(SVector{D,Int})
    else
        setindex(exp, new_exp_val, idx)
    end
    
    idxs_tail = tail(idxs)
    if length(idxs_tail) > 0
        return derivative(new_exp,val/h*exp_val,h,idxs_tail...)
    else
        return val/h*exp_val, new_exp
    end
end
"""
    derivative(m1::Monomial{T,D}, 
        h::T, idxs::Vararg{Int,N}) where {T,D,N} -> Monomial

Compute the derivative of a monomial with respect to a given index.
The monomials can be scaled as `((x-x0)/h)^n`, so 
the value h needs to be passed to compute the derivative.

# Arguments
- `m1::Monomial{T,D}`: Monomial to compute the derivative.
- `h::T`: Value of h.
- `idxs::Vararg{Int,N}`: Indexes to compute the derivative.

# Example 
```julia
m1 = Monomial(2.0, SVector(1,0))  # 2x

h = 0.5 
∂(m1, h, 1)  # 2/h
"""
function derivative(m1::Monomial{T,D}, h::S, idxs::Vararg{Int,N}) where {T,D,S,N}
    val,exp =  derivative(m1.exp, m1.val, h, idxs...)
    Monomial(val,exp)
end

function ∇(m1::Monomial{T,D},h::S = 1.0) where {T,D,S} 
    grad = SVector{D}(derivative(m1,h,i) for i in 1:D)
    return grad
end

function ∇(sm1::AbstractVector{Monomial{T,D}},h::S = 1.0) where {T,D,S} 
    U = length(typeof(sm1)) #does only work for statically sized vectors
    # grad = Mat{U,D}((i,j)->derivative(sm1[i],h,j))
    grad = SMatrix{U,D}(derivative(sm1[i],h,j) for i in 1:U, j in 1:D)
    return grad
end

const ∂ = derivative



Base.:*(m1::Monomial{T,D}, m2::Monomial{T,D}) where {T,D} = 
    Monomial(m1.val * m2.val, m1.exp + m2.exp)

Base.:*(a::Number, m::Monomial{T,D}) where {T,D} = 
    Monomial{T,D}(a * m.val, m.exp)
Base.:*(m::Monomial{T,D}, a::Number) where {T,D} = 
    Monomial{T,D}(m.val * a, m.exp)

Base.:-(m::Monomial) = Monomial(-m.val, m.exp)
Base.:+(m::Monomial) = m

Base.:/(m::Monomial{T,D}, a::Number) where {T,D} = 
    Monomial{T,D}(m.val / a, m.exp)


# function LinearAlgebra.dot(v::AbstractVector{Monomial{T,D}},n::AbstractVector{T}) where {T,D}
#     return sum(v[i]*n[i] for i in 1:D)
# end


########## Display Functionality ##############

# Simple string representation for monomials
function _monomial_string(m::Monomial{T,D}) where {T,D}
    # Handle zero case
    if m.val == 0
        return "0"
    end
    
    # Handle constant term
    if all(m.exp .== 0)
        return string(round(m.val, sigdigits=3))
    end
    
    # Build the string
    result = ""
    
    # Add coefficient (skip if 1, show minus if -1)
    if m.val == 1
        # Don't add anything for coefficient 1
    elseif m.val == -1
        result = "-"
    else
        result = string(round(m.val, sigdigits=3))
    end
    
    # Add variables with exponents
    for (i, exp) in enumerate(m.exp)
        if exp > 0
            if exp == 1
                result *= "x" * string(i)
            else
                result *= "x" * string(i) * "^" * string(exp)
            end
        end
    end
    
    return result
end

# Display methods
function Base.show(io::IO, m::Monomial{T,D}) where {T,D}
    print(io, _monomial_string(m))
end

function Base.show(io::IO, ::MIME"text/plain", m::Monomial{T,D}) where {T,D}
    println(io, _monomial_string(m))
end

"""
    is_constant(m::Monomial{T,D}) where {T,D}

Check if a monomial is constant (all exponents are zero).

# Examples
```julia
m1 = Monomial(5.0, SVector(0, 0))  # constant
m2 = Monomial(2.0, SVector(1, 0))  # not constant
is_constant(m1)  # true
is_constant(m2)  # false
```
"""
function is_constant(m::Monomial{T,D}) where {T,D}
    return all(m.exp .== 0)
end


