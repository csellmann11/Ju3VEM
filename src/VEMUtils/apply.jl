#################
# Boundary conditions
#################
function mean_diag(A::SparseMatrixCSC{T}) where T<:Real
    @assert A.m == A.n "n_rows must be n_cols"
    s = zero(T)
    for i in axes(A,1)
        s += A[i,i]
    end

    return s/A.m
end

struct RhsData{T<:Float64}
    rhs_vec::Vector{T}
    md     ::T
end


"""
    apply!(A, f, ch; get_rhs_data=false)

Apply Dirichlet and Neumann boundary conditions to the system matrix `A`
and the right-hand side vector `f` using the constraint handler `ch`.

# Arguments
- `A::SparseMatrixCSC{T}`: A square sparse matrix.
- `f::AbstractVector{T}`: The right-hand side vector.
- `ch::ConstraintHandler{D, U}`: Contains Dirichlet (`d_bcs`) and Neumann (`n_bcs`) BCs.
- `get_rhs_data::Bool`: Optional flag. If `true`, collects modifications made to `f`.

# Returns
- An instance of `RhsData{T}` containing:
  - `rhs_vec`: The vector of right-hand side modifications.
  - `md`: The mean of the diagonal elements of `A`.
"""
function apply!(A::SparseMatrixCSC{T}, f::AbstractVector{T},
                ch::ConstraintHandler{D, U};zero_dbcs = false,
                get_rhs_data::Bool = false) where {T<:Float64, D, U}

    @assert A.m == A.n "Matrix A must be square."

    # Compute the mean of the diagonal elements.
    mean_diag_val = mean_diag(A)
    
    # Initialize the right-hand side data.
    rhs_vector = get_rhs_data ? zeros(T, size(f)) : T[]
    rhs_data = RhsData{T}(rhs_vector, mean_diag_val)
    
    # Extract Dirichlet BCs and precompute a Boolean mask.
    dirichlet_bcs = ch.d_bcs
    n = size(A, 1)
    is_dirichlet = falses(n)
    for key in keys(dirichlet_bcs)
        if key <= n
            is_dirichlet[key] = true
        end
    end

    # Retrieve sparse matrix data.
    row_indices = rowvals(A)
    values = nonzeros(A)
    
    # Loop over each column to apply Dirichlet conditions.
    for col in axes(A, 2)
        col_has_bc = is_dirichlet[col]
        # Only perform dictionary lookup if the column has a BC.
        constraint_value = col_has_bc ? (dirichlet_bcs[col]* !zero_dbcs) : nothing

        for i in nzrange(A, col)
            row = row_indices[i]
            current_val = values[i]

            if is_dirichlet[row] && row == col
                # Diagonal entry: set to mean_diag_val and adjust f.
                values[i] = mean_diag_val
                f[row] = constraint_value * mean_diag_val
            elseif is_dirichlet[row]
                # Off-diagonal entry in a Dirichlet row: zero it out.
                values[i] = zero(T)
            elseif col_has_bc
                # Non-Dirichlet row influenced by a Dirichlet column.
                f[row] -= constraint_value * current_val
                if get_rhs_data
                    rhs_data.rhs_vec[row] = constraint_value * current_val
                end
                values[i] = zero(T)
            end
        end
    end

    # Apply Neumann boundary conditions.
    for (idx, neumann_value) in ch.n_bcs
        if idx <= n
            f[idx] += neumann_value
        end
    end

    return rhs_data
end




function apply_rhs!(rhs::Vector{Float64},
    ch::ConstraintHandler{D,U},rhs_data::RhsData) where {D,U}


    rhs .-= rhs_data.rhs_vec

    md = rhs_data.md
    for (idx,val) in ch.d_bcs
        rhs[idx] = val*md
    end

    for (idx,val) in ch.n_bcs
        rhs[idx] += val
    end
end