
using Base:@propagate_inbounds
import Base: getindex, setindex!, size, similar, convert
"""
Type providing a re-sizable array.

The size of a `CachedArray` is changed via the [`setsize!`](@ref) function.

A `CachedArray` can be build with the constructors
- [`CachedArray(a::AbstractArray)`](@ref)
- [`CachedArray(T,N)`](@ref)

```jldoctests
using Gridap.Arrays
# Create an empty CachedArray
a = CachedArray(Float64,2)
# Resize to new shape (2,3)
setsize!(a,(2,3))
size(a)
# output
(2, 3)
```
"""
mutable struct CachedArray{T,N,A<:AbstractArray{T,N}} <: AbstractArray{T,N}
  array::A
  const buffer::Dict{NTuple{N,Int},A}

  @doc """
      CachedArray(a::AbstractArray)
  
  Constructs a `CachedArray` from a given array.
  """
  function CachedArray(array::A) where {T,N,A<:AbstractArray{T,N}}
    buffer = Dict{NTuple{N,Int},A}()
    buffer[size(array)] = array
    new{T,N,A}(array,buffer)
  end
end

function CachedArray{T,N,A}() where {T,N,A<:AbstractArray{T,N}}
  array = Array{T,N}(undef,ntuple(i->0,Val(N))...)
  CachedArray(array)
end

"""
    const CachedMatrix{T,A} = CachedArray{T,2,A}
"""
const CachedMatrix{T,A} = CachedArray{T,2,A}
const DefCachedMatrix{T} = CachedArray{T,2,Matrix{T}}

"""
    const CachedVector{T,A} = CachedArray{T,1,A}
"""
const CachedVector{T,A} = CachedArray{T,1,A}
const DefCachedVector{T} = CachedVector{T,Vector{T}}


CachedVector(a::AbstractVector) = CachedArray(a)


CachedMatrix(a::AbstractMatrix) = CachedArray(a)

"""
    CachedArray(T,N)

Constructs an empty `CachedArray` of element type `T` and `N` dimensions.
"""
function CachedArray(T,N)
  s = tuple([0 for i in 1:N]...)
  a = Array{T,N}(undef,s)
  CachedArray(a)
end


function CachedVector(T)
  CachedArray(T,1)
end


function CachedMatrix(T)
  CachedArray(T,2)
end

size(self::CachedArray) = size(self.array)

"""
Changes the size of the `CachedArray` `a` to the size described the the tuple
`s`.
After calling `setsize!`, the array can store uninitialized values.
"""
function setsize!(a::CachedArray{T,N},s::NTuple{N,Int},
  set_zero::Bool = true, init_fun::Function = identity) where {T,N}

  if s != size(a.array) 
    a.array = get!(a.buffer,s) do
      A = zeros(T,s)
      init_fun(A)
    end
  end

  if set_zero a.array .= Ref(zero(T)) end
  return nothing
end


function setsize!(a::CachedArray{T,N},s::NTuple{N,<:Integer}) where {T,N}
  _s::NTuple{N,Int} = s
  setsize!(a,_s)
end

function setsize!(v::CachedVector{T},s::Integer) where {T}
  setsize!(v,(s,))
end


@propagate_inbounds function getindex(self::CachedArray, kj::Integer)
    self.array[kj]
end

@propagate_inbounds function getindex(self::CachedArray{T,N}, kj::Vararg{Integer,N}) where {T,N}
    self.array[kj...]
end

@propagate_inbounds function setindex!(B::CachedArray, v, kj::Integer)
    B.array[kj] = v
    v
end

@propagate_inbounds function setindex!(B::CachedArray{T,N}, v, kj::Vararg{Integer,N}) where {T,N}
    B.array[kj...] = v
    v
end

function similar(::Type{CachedArray{T,N,A}},s::Tuple{Vararg{Int}}) where {T,N,A}
  a = similar(A,s)
  CachedArray(a)
end

Base.convert(::Type{CachedArray{T,N,A}},a::CachedArray{T,N,A}) where {T,N,A} = a

function Base.convert(::Type{CachedArray{T,N,A}},a::CachedArray) where {T,N,A}
  array = convert(A,a.array)
  CachedArray(array)
end