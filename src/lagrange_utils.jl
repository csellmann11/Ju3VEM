function lagrange_shape_function(idx::Int, coord::Real, ::Val{O}, ::Val{U}) where {O, U}
    # Check if coord is in the valid range
    if !(coord >= -1 && coord <= 1)
        error("Coordinate must be in the range [-1, 1], got $coord")
    end
    
    # For U=1, return the scalar directly
    if U == 1
        return _get_scalar_basis(idx, coord, Val(O))
    end
    
    # For U>1, determine which component should be non-zero
    component = (idx - 1) % U + 1
    basis_idx = (idx - 1) ÷ U + 1
    
    # Create the vector with the appropriate component set
    result = zeros(SVector{U, Float64})
    
    # Get the scalar basis function for the correct index
    scalar_val = _get_scalar_basis(basis_idx, coord, Val(O))
    
    # Set the appropriate component
    result = setindex(result, scalar_val, component)
    
    return result
end


function _get_scalar_basis(idx::Int, coord::Real, ::Val{0})
    if idx == 1
        return 1.0
    else
        error("Invalid index for first order Lagrange shape function. Must be 1 or 2.")
    end
end

# Helper function to compute scalar basis functions
function _get_scalar_basis(idx::Int, coord::Real, ::Val{1})
    if idx == 1
        return 0.5 * (1 - coord)
    elseif idx == 2
        return 0.5 * (1 + coord)
    else
        error("Invalid index for first order Lagrange shape function. Must be 1 or 2.")
    end
end

function _get_scalar_basis(idx::Int, coord::Real, ::Val{2})
    if idx == 1
        return 0.5 * coord * (coord - 1)
    elseif idx == 2
        return 1 - coord^2
    elseif idx == 3
        return 0.5 * coord * (coord + 1)
    else
        error("Invalid index for second order Lagrange shape function. Must be 1, 2, or 3.")
    end
end

function _get_scalar_basis(idx::Int, coord::Real, ::Val{3})
    if idx == 1
        return -1/16 * (coord - 1) * (9*coord^2 - 1)
    elseif idx == 2
        return 9/16 * (coord + 1) * (1 - 3*coord) * (1 - coord)
    elseif idx == 3
        return 9/16 * (coord + 1) * (1 + 3*coord) * (coord - 1)
    elseif idx == 4
        return 1/16 * (coord + 1) * (9*coord^2 - 1)
    else
        error("Invalid index for third order Lagrange shape function. Must be between 1 and 4.")
    end
end

function _get_scalar_basis(idx::Int, coord::Real, ::Val{4})
    if idx == 1
        return 1/6 * (coord - 1) * coord * (coord - 1/2) * (coord + 1/2)
    elseif idx == 2
        return -4/3 * (coord + 1) * coord * (coord - 1/2) * (coord - 1)
    elseif idx == 3
        return 1 - 5/2 * coord^2 + 3/2 * coord^4
    elseif idx == 4
        return -4/3 * (coord + 1) * coord * (coord + 1/2) * (coord - 1)
    elseif idx == 5
        return 1/6 * (coord + 1) * coord * (coord + 1/2) * (coord - 1/2)
    else
        error("Invalid index for fourth order Lagrange shape function. Must be between 1 and 5.")
    end
end

# Convenience function for scalar output
function lagrange_shape_function(idx::Int, coord::Real, ::Val{O}) where {O}
    return lagrange_shape_function(idx, coord, Val(O), Val(1))
end


function interpolate_edge(ξ::T,n1::V,n2::V) where {T,V<:AbstractVector{<:Real}}
    return _get_scalar_basis(1,ξ,Val(1))*n1 + _get_scalar_basis(2,ξ,Val(1))*n2
end