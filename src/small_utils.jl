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



 
function get_unique_values(v::AbstractVector{<:AbstractVector{Int}},
    ids::AbstractVector{Int})

    @assert maximum(ids) <= length(v) "ids are out of bounds"
    max_len = sum(length,v[i] for i in ids)
    pointer = 0
    unique_values = @no_escape begin 
        _unique_values = @alloc(Int,max_len)
        @inbounds for id in ids
            values = v[id]
            for value in values
                if value ∉ @view(_unique_values[1:pointer])
                    _unique_values[pointer+1] = value
                    pointer += 1
                end
            end
        end
        _unique_values[1:pointer]
    end
    return unique_values
end

"""
    get_unique_values(v::AbstractVector{<:AbstractVector{Int}})
# Description:
-   returns the unique ids in the vector
# Arguments:
`v::AbstractVector{<:AbstractVector{Int}}`: the vector to get the unique ids from
# Returns:
 `Vector{Int}`: the unique ids in the vector
# Example:
    get_unique_values([[1,2,3],[4,5,6],[1,2,3]]) == [1,2,3,4,5,6]
"""
function get_unique_values(v::AbstractVector{<:AbstractVector{Int}})

    return get_unique_values(v,1:length(v))
end


function find_single_intersec(vecs_col::AbstractVector{<:AbstractVector{Int}},
    ids::AbstractVector{Int})
    for a in vecs_col[first(ids)]
        all(a in vecs_col[i] for i in @views ids[2:end]) && return a
    end
    return 0

end
"""
    find_single_intersec(vecs::AbstractVector{<:AbstractVector{Int}})
# Description:
-   returns the single intersection of the vectors
-   returns the first intersection found if there are multiple
-   returns 0 if no intersection is found
# Arguments:
`vecs::AbstractVector{<:AbstractVector{Int}}`: the vectors to find the single intersection of
# Returns:
`Int`: the single intersection of the vectors
# Example:
    find_single_intersec([[1,2,3],[4,5,6],[1,2,3]]) == 1
"""
function find_single_intersec(vecs::AbstractVector{<:AbstractVector{Int}}) 
    return find_single_intersec(vecs,1:length(vecs))

end
