using StaticArrays
using Chairmarks
using Test
include("monomials.jl")
include("polynomials.jl")
include("../stretched_matrices.jl")


bi = BaseInfo{2,2,1}()
base = get_base(bi)
base
p = Polynomial(ones(length(base)),base.base)

x = SA[1.0,2.0]
p(x)
expected_val = 1 + x[1] + x[2] + x[1]^2 + x[2]^2 + x[1]*x[2]
@test p(x) ≈ expected_val



expected_grad = [SA[1.0,2.0,1.0],SA[1.0,1.0,2.0]]
@test ∇p(p,1.0)[1].coeffs ≈ expected_grad[1]
@test ∇p(p,1.0)[2].coeffs ≈ expected_grad[2]


exptected_eval = SA[1.0 + 2.0*x[1] + 1.0*x[2], 1.0 + 1.0*x[1] + 2.0*x[2]]
@test ∇p(p,1.0)(x) ≈ exptected_eval


pmat = SA[p p; p p]

@test ∇p(pmat,1.0)[1,1,1] == ∇p(p,1.0)[1] 

@test ∇x(p,1.0,x) ≈ ∇p(p,1.0)(x)
