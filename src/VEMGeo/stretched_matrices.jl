
struct StretchedMatrix{T,λ,M<:AbstractMatrix{T}} <: AbstractMatrix{T}
    data::M
end
function StretchedMatrix{λ}(A::AbstractMatrix{T}) where {T,λ}
    return StretchedMatrix{T,λ,typeof(A)}(A)
end

# Hilfsfunktion zum Erstellen einer StretchedMatrix
stretch(A::AbstractMatrix, ::Val{λ} = Val(1)) where λ = StretchedMatrix{λ}(A)

# Optimized getindex-Methode für StretchedMatrix
Base.@propagate_inbounds function Base.getindex(A::StretchedMatrix{T,λ,M}, i::Int, j::Int) where {T,λ,M}   
    # This condition is already efficient and clear.
    if (i - 1) % λ == (j - 1) % λ
        # Use integer division, which is generally faster than floating-point division + ceil.
        orig_i = (i - 1) ÷ λ + 1
        orig_j = (j - 1) ÷ λ + 1
        return A.data[orig_i, orig_j]
    else
        return zero(T)
    end
end

# Optimized getindex-Methode für StretchedMatrix
Base.@propagate_inbounds function Base.getindex(A::StretchedMatrix{T,1,M}, i::Int, j::Int) where {T,M}   
    return A.data[i,j]
end

Base.size(A::StretchedMatrix{T,λ}) where {T,λ} = size(A.data) .* λ

function Base.adjoint(A::StretchedMatrix{T,λ}) where {T,λ}
    return StretchedMatrix{λ}(A.data')
end

function Base.copy(A::StretchedMatrix{T,λ}) where {T,λ}
    StretchedMatrix{T,λ,typeof(A.data)}(copy(A.data))
end

function destretch!(B,BS::StretchedMatrix{Float64,N}) where {N}
    @turbo B .= 0.0
    for idx in CartesianIndices(BS.data)
        i,j = idx.I 
        @turbo for k in 1:N
            idx1 = (i-1)*N + k;
            idx2 = (j-1)*N + k 
            B[idx1,idx2] = BS.data[idx]
        end
    end
    B
end

function destretch!(B,BS::StretchedMatrix{Float64,1})
    @turbo B .= BS.data
    B
end

function destretch(BS::StretchedMatrix{Float64,N}) where {N}
    B = zeros(Float64, size(BS)...)
    destretch!(B,BS)
end


@inline destretch(BS) = BS
function generic_spaced_matmul!(C::AbstractMatrix{T},
    AS::AbstractMatrix{T},BS::AbstractMatrix{T}) where {T<:Real}

    A = destretch(AS); B = destretch(BS)
    Octavian.matmul!(C,A,B)
    return nothing
end