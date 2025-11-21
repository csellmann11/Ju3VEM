"""
    get_next_idx(v::AbstractVector,i::Int)
# Description:
-   returns the index of the next element in the vector
# Arguments:
`v::AbstractVector`: the vector to get the next index from
`i::Int`: the index to get the next index from
# Returns:
`Int`: the index of the next element in the vector
# Example:
    get_next_idx([1,2,3],1) == 2
    get_next_idx([1,2,3],2) == 3
    get_next_idx([1,2,3],3) == 1
"""
function get_next_idx(v::AbstractVector,i::Int)
    return i == length(v) ? 1 : i + 1
end

"""
    get_prev_idx(v::AbstractVector,i::Int)
# Description:
-   returns the index of the previous element in the vector
# Arguments:
`v::AbstractVector`: the vector to get the previous index from
`i::Int`: the index to get the previous index from
# Returns:
`Int`: the index of the previous element in the vector
"""
function get_prev_idx(v::AbstractVector,i::Int)
    return i == 1 ? length(v) : i - 1
end



 
function get_unique_values(v::AbstractVector{<:AbstractVector{I}},
    ids::AbstractVector{<:Integer}) where {I<:Integer}

    @assert maximum(ids) <= maximum(keys(v)) "ids are out of bounds"
    max_len = sum(length,v[i] for i in ids)
    pointer = 0
    _unique_values = ID_VEC_TYPE{I}(undef, max_len)
    for id in ids
        values = v[id]
        for value in values
            if value ∉ @view(_unique_values[1:pointer])
                _unique_values[pointer+1] = value
                pointer += 1
            end
        end
    end
    return _unique_values[1:pointer]
    
    
end

"""
    get_unique_values(v::AbstractVector{<:AbstractVector{<:Integer}})
# Description:
-   returns the unique ids in the vector
# Arguments:
`v::AbstractVector{<:AbstractVector{<:Integer}}`: the vector to get the unique ids from
# Returns:
 `Vector{<:Integer}`: the unique ids in the vector
# Example:
    get_unique_values([[1,2,3],[4,5,6],[1,2,3]]) == [1,2,3,4,5,6]
"""
function get_unique_values(v::AbstractVector{<:AbstractVector{<:Integer}})

    return get_unique_values(v,1:length(v))
end


function find_single_intersec(vecs::Vararg{AbstractVector{<:Integer},N}; 
    ids::AbstractVector{<:Integer} = 1:N) where N
    for a in vecs[first(ids)]
        all(a in vecs[i] for i in @views ids[2:end]) && return a
    end
    return 0

end


function find_single_intersec(vecs_col::AbstractVector{<:AbstractVector{<:Integer}},
    ids::AbstractVector{<:Integer} = 1:length(vecs_col))
    for a in vecs_col[first(ids)]
        all(a in vecs_col[i] for i in @views ids[2:end]) && return a
    end
    return 0

end
"""
    find_single_intersec(vecs::AbstractVector{<:AbstractVector{<:Integer}})
# Description:
-   returns the single intersection of the vectors
-   returns the first intersection found if there are multiple
-   returns 0 if no intersection is found
# Arguments:
`vecs::AbstractVector{<:AbstractVector{<:Integer}}`: the vectors to find the single intersection of
# Returns:
`Int`: the single intersection of the vectors
# Example:
    find_single_intersec([[1,2,3],[4,5,6],[1,2,3]]) == 1
"""
# function find_single_intersec(vecs::AbstractVector{<:AbstractVector{Int}}) 
#     return find_single_intersec(vecs,1:length(vecs))

# end



function max_node_distance(nodes::AbstractVector{V},
    ids::AbstractVector{<:Integer}) where {T<:Real,V <: AbstractVector{T}}
    max_dist = zero(T)
    n_nodes = length(ids) 
    @inbounds for i in 1:n_nodes-1, j in i+1:n_nodes
        n1 = nodes[ids[i]]; n2 = nodes[ids[j]]
        dist = norm(n2-n1)
        max_dist = max(max_dist, dist)
    end
    return max_dist
end



function max_node_distance(nodes::AbstractVector) 
    return max_node_distance(nodes,1:length(nodes))
end



function resize_with_trailing_zeros!(v::AbstractVector{T},size::Integer) where T
    current_size = length(v)
    resize!(v,size)
    v[current_size+1:end] .= zero(T)
    v
end

function resize_and_fill!(dest::AbstractVector,data)
    resize!(dest,size(data)|>only)
    for i in eachindex(dest,data)
        dest[i] = data[i]
    end
end

function resize_and_fill!(dest::AbstractVector,
    data,ids::AbstractVector{<:Integer})
    resize!(dest,size(ids)|>only)
    for i in eachindex(dest,ids)
        dest[i] = data[ids[i]]
    end
end


function get_unscaled_normal(n1::StaticVector{2,Float64},n2::StaticVector{2,Float64})
    return @inbounds SVector(n1[2]-n2[2],n2[1]-n1[1])
end

"""
function returns the scaled normal. 
The second node is passed first for a ccw ordering
"""
function get_normal(n1::V,n2::V) where {T<:Number,V <:AbstractVector{T}}
    n = get_unscaled_normal(n1,n2)
    len = sqrt(sum(abs2,n))
    return n/len, len
end
