using Symbolics
using Symbolics.RuntimeGeneratedFunctions: @RuntimeGeneratedFunction


"""
    build_ifelse_chain(conditions, actions, else_action = nothing)

Build an if-else chain from a list of conditions and actions.

# Arguments
- `conditions::Vector{Expr}`: List of conditions to check.
- `actions::Vector{Expr}`: List of actions to execute if the condition is true.
- `else_action::Expr`: Action to execute if no condition is true.

# Returns
- `Expr`: An if-else chain expression.
"""
function _build_ifelse_chain(conditions, actions, else_action = nothing)

    @assert length(conditions) == length(actions)

    if else_action === nothing
        else_action = Expr(:call, :throw, Expr(:call, :ErrorException, "No matching condition found"))
    end

    # Start with the final else
    expr = Expr(:block, else_action)
    
    # Work backwards through conditions
    for i in reverse(eachindex(conditions))
        expr = Expr(:if, conditions[i], Expr(:block, actions[i]), expr)
    end
    
    return expr
end

"""
    _generate_condition_matrix(rows, cols)

Generate a matrix of conditions for an if-else chain.
Matrix at index (a,b) has the condition: i == (a-1) && j == (b-1)

# Arguments
- `rows::Int`: Number of rows in the matrix.
- `cols::Int`: Number of columns in the matrix.
"""
function _generate_condition_matrix(dims::Vararg{Int,N}) where N
    conditions = Array{Expr}(undef, dims...)

    for idxs in CartesianIndices(dims)
        _idxs = idxs.I 
        # expr = :(_i1 == $(_idxs[1]-1) && _i2 == $(_idxs[2]-1) && _i3 == $(_idxs[3]-1))
        expr = Expr(:(&&))
        for i in 1:N
            name = Symbol(:_i,i)
            val  = _idxs[i] - 1
            push!(expr.args, Expr(:call, :(==), name, val))
        end
        conditions[idxs] = expr
    end
    return conditions
end


"""
    _generate_polyprod_expr_matrix(n::Int)

Generate a matrix of expressions for the product of two polynomials.
Matrix at index (a,b) has the expression: p_i * p_j

# Arguments
- `n::Int`: Order of the polynomials.
"""
function _generate_polyprod_expr_matrix2d_to_3d(n::Int) 
    base    = get_base(BaseInfo{3,1,1}()).base 
    coeffs1 = Symbolics.scalarize(Symbolics.variables(:c,1:4)) #|> SVector{4,Num}
    coeffs2 = Symbolics.scalarize(Symbolics.variables(:d,1:4)) #|> SVector{4,Num}
    coeffs  = vcat(coeffs1,coeffs2)

    lin1_poly = Polynomial(coeffs1,base);
    lin2_poly = Polynomial(coeffs2,base);

    expr_mat = Matrix{Expr}(undef,n,n)
    for i in axes(expr_mat,1)
        p_i = poly_pow(lin1_poly,Val(i-1))
        for j in axes(expr_mat,2)
            p_j = poly_pow(lin2_poly,Val(j-1))
            p_prod = p_i * p_j           
            f_expr = Symbolics.build_function(p_prod.coeffs,coeffs,boundscheck= true, cse = true)[2]
            expr_mat[i,j] = f_expr.args[2]
        end
    end
    return expr_mat
end


"""
    _generate_polyprod_expr_matrix2d_to_2(n::Int)

Prepares the coefficients for (x + const1)^i * (x + const2)^j
"""
function _generate_polyprod_expr_matrix2d_to_2(n::Int) 
    base    = get_base(BaseInfo{2,1,1}()).base 
    coeff1 = Symbolics.scalarize(Symbolics.variable(:c,1)) #|> SVector{4,Num}
    coeff2 = Symbolics.scalarize(Symbolics.variable(:d,1)) #|> SVector{4,Num}
    coeffs  = [coeff1,coeff2]

    lin1_poly = Polynomial([coeff1,Num(1),Num(0)],base)
    lin2_poly = Polynomial([coeff2,Num(0),Num(1)],base)


    expr_mat = Matrix{Expr}(undef,n,n)
    for i in axes(expr_mat,1)
        p_i = poly_pow(lin1_poly,Val(i-1))
        for j in axes(expr_mat,2)
            p_j = poly_pow(lin2_poly,Val(j-1))
            p_prod = p_i * p_j           
            f_expr = Symbolics.build_function(p_prod.coeffs,coeffs,boundscheck= true, cse = true)[2]
            expr_mat[i,j] = f_expr.args[2]
        end
    end
    return expr_mat
end


function _generate_polyprod_expr_matrix3d_to_2d(n::Int) 
    base    = get_base(BaseInfo{2,1,1}()).base 
    coeffs1 = Symbolics.scalarize(Symbolics.variables(:c,1:3)) #|> SVector{3,Num}
    coeffs2 = Symbolics.scalarize(Symbolics.variables(:d,1:3)) #|> SVector{3,Num}
    coeffs3 = Symbolics.scalarize(Symbolics.variables(:e,1:3)) #|> SVector{3,Num}
    coeffs  = vcat(coeffs1,coeffs2,coeffs3)


    lin1_poly = Polynomial(coeffs1,base);
    lin2_poly = Polynomial(coeffs2,base);
    lin3_poly = Polynomial(coeffs3,base);

    expr_mat = Array{Expr}(undef,n,n,n)
    for i in axes(expr_mat,1)
        p_i = poly_pow(lin1_poly,Val(i-1))
        for j in axes(expr_mat,2)
            p_j = poly_pow(lin2_poly,Val(j-1))
            for k in axes(expr_mat,3)
                p_k = poly_pow(lin3_poly,Val(k-1))
                p_prod = p_i * p_j * p_k           
                f_expr = Symbolics.build_function(p_prod.coeffs,coeffs,checkbounds = true, cse = true)[2]
                expr_mat[i,j,k] = f_expr.args[2]
            end
        end
    end
    return expr_mat
end

@generated function compute_transformation_coeffs2d_to_3d!(ˍ₋out::AbstractVector{T},
    m::Monomial{T,2},
    bc::StaticVector{2},
    ::Val{N} = Val(3)) where {T,N}

    expr_mat = _generate_polyprod_expr_matrix2d_to_2(N+1)
    cond_mat = _generate_condition_matrix(size(expr_mat)...)
    if_else_expr = _build_ifelse_chain(cond_mat, expr_mat)

    return quote 
        ˍ₋arg1 = bc 
        _i1,_i2 = m.exp
        $if_else_expr

        ˍ₋out .*= m.val
        return ˍ₋out
    end

end

function compute_transformation_coeffs2d_to_2d(
    m::Monomial{T,2},bc::StaticVector{2},::Val{N} = Val(3)) where {T,N}
    len = get_base_len(2,sum(m.exp),1)
    
    ˍ₋out = FixedSizeVector{Float64}(undef,len)
    compute_transformation_coeffs2d_to_2d!(ˍ₋out,m,bc,Val(N))
    return ˍ₋out
end

@generated function compute_transformation_coeffs3d_to_2d!(ˍ₋out::AbstractVector{T},
    m::Monomial{T,3},p0::StaticVector,u::StaticVector,v::StaticVector,::Val{N} = Val(3)) where {T,N}

    expr_mat = _generate_polyprod_expr_matrix3d_to_2d(N+1)
    cond_mat = _generate_condition_matrix(size(expr_mat)...)
    if_else_expr = _build_ifelse_chain(cond_mat, expr_mat)

    return quote 
        ˍ₋arg1      = hcat(p0,u,v)'[:]
        _i1,_i2,_i3 = m.exp
        $if_else_expr

        ˍ₋out       .*= m.val
        return ˍ₋out
    end

end

function compute_transformation_coeffs3d_to_2d(
    m::Monomial{T,3},p0::StaticVector{3},u::StaticVector,v::StaticVector,::Val{N} = Val(3)) where {T,N}

    order       = sum(m.exp)
    len         = get_base_len(2,order,1)
    ˍ₋out       = FixedSizeVector{Float64}(undef,len) 
    compute_transformation_coeffs3d_to_2d!(ˍ₋out,m,p0,u,v,Val(N)) 
    return ˍ₋out 
end


# #TODO: add safety check for polynomial order, reduce number of computed entries in expression matrix
# @generated function compute_transformation_coeffs3d_to_2d(
#     m::Monomial{T,3},u::StaticVector,v::StaticVector,::Val{N} = Val(3)) where {T,N}
#     expr_mat = _generate_polyprod_expr_matrix3d_to_2d(N)
#     cond_mat = _generate_condition_matrix(size(expr_mat)...)
#     if_else_expr = _build_ifelse_chain(cond_mat, expr_mat)

#     return quote 
#         ˍ₋arg1      = hcat(u,v)'[:]
#         order       = sum(m.exp)
#         len         = get_base_len(2,order,1)
#         ˍ₋out       = FixedSizeVector{Float64}(undef,len) 
#         _i1,_i2,_i3 = m.exp
#         $if_else_expr

#         ˍ₋out       .*= m.val
#         return ˍ₋out
#     end
# end


#TODO: add safety check for polynomial order, reduce number of computed entries in expression matrix
@generated function compute_transformation_coeffs2d_to_3d(
    m::Monomial{T,2},p0::StaticVector{2},u::StaticVector,v::StaticVector,::Val{N} = Val(3)) where {T,N}
    expr_mat = _generate_polyprod_expr_matrix2d_to_3d(N)
    cond_mat = _generate_condition_matrix(size(expr_mat)...)
    if_else_expr = _build_ifelse_chain(cond_mat, expr_mat)

    return quote 
        ˍ₋arg1      = vcat(SA[p0[1]],u,SA[p0[2]],v)
        order       = sum(m.exp) 
        len         = get_base_len(3,order,1)
        ˍ₋out       = FixedSizeVector{Float64}(undef,len) 
        _i1,_i2    = m.exp
        $if_else_expr

        ˍ₋out       .*= m.val
        return ˍ₋out
    end
end







