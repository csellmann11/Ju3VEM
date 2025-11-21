###########################################################################
###########################################################################
# FlattenVecs
###########################################################################
###########################################################################

""" 
    mutable struct FlattenVecs{N,T,V<:AbstractVector{T}} <: AbstractVector{T}
        v::ApplyArray{T, 1, typeof(vcat), NTuple{N,V}}
    end
FlattenVecs stores a statically known number of vectors. 
Each vector can be resizes later on. 

PERFORMANCE WARNING: Indexing is slower than for standard arrays

Indexing: 
```julia 
    a = [1,2]
    b = [3,4] 
    v = FlattenVecs(a,b)

    # Indexing:
    v[3] = 3 
    v.v[1][1] = 3

    push!(a,5)

    length(v) # returns 5
    
```

"""
mutable struct FlattenVecs{N,T,V<:AbstractVector{T}} <: AbstractVector{T}
    # uses LazyArrays under the hood
    v::ApplyArray{T, 1, typeof(vcat), NTuple{N,V}}
end

function FlattenVecs(v::Vararg{<:AbstractVector{T},N}) where {N,T}
    return FlattenVecs(Vcat(v...))
end

@inline FlattenVecs{N,T,V}() where {N,T,V<:AbstractVector{T}} = 
        FlattenVecs(ntuple(i -> V(), Val(N))...)

@inline FlattenVecs{N,T}() where {N,T} = 
    FlattenVecs{N,T,Vector{T}}()


Base.@propagate_inbounds Base.getindex(t::FlattenVecs, idx::Int) = t.v[idx]  
Base.@propagate_inbounds Base.setindex!(t::FlattenVecs, val, idx::Int) = t.v[idx] = val

Base.length(t::FV) where FV<:FlattenVecs = length(t.v)
Base.size(t::FV) where FV<:FlattenVecs = (length(t),)

Base.eltype(::FlattenVecs{N,T})  where {N,T} = T


@inline Base.@propagate_inbounds get_first(v::FlattenVecs)  = v.v.args[1]
@inline Base.@propagate_inbounds get_second(v::FlattenVecs) = v.v.args[2]
@inline Base.@propagate_inbounds get_third(v::FlattenVecs)  = v.v.args[3]
@inline Base.@propagate_inbounds get_fourth(v::FlattenVecs) = v.v.args[4]


