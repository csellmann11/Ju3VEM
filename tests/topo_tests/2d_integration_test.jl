using StaticArrays, WriteVTK
using OrderedCollections, Bumper
using LinearAlgebra, Statistics
using SmallCollections, Chairmarks
using LoopVectorization
using Test

# Include the main module (professional approach)
include("../../src/small_utils.jl")
include("../../src/topo.jl")
include("../../src/vtkexports.jl")
include("../../src/element_refinement.jl")
include("../../src/lagrange_utils.jl")
include("../../src/element_coarsening.jl")
include("../../src/triangulation.jl")

include("../../src/Polynomials/polynomials.jl")
include("../../src/Polynomials/monomials.jl")
include("../../src/polynomial_integration.jl")

@testset "2D Integration Tests" begin
    node_coords = 2 .* [SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0), SVector(0.0, 1.0)]

    m1 = Monomial(1.0, SVector(0, 0))

    expected_area = 4.0
    @test isapprox(integral(m1, node_coords), expected_area; atol=1e-10)


    m2 = Monomial(1.0, SVector(1, 0))
    bc = SVector(1.0, 1.0)
    integral(m2, node_coords, 1.0, bc)

    expected_val = 0.0
    @test isapprox(integral(m2, node_coords, 1.0, SVector(1.0, 1.0)), expected_val; atol=1e-10)


    expected_val = bc[1]
    @test isapprox(integral(m2, node_coords, 1.0) / expected_area, expected_val; atol=1e-10)

    m3 = Monomial(1.0, SVector(0, 1))
    expected_val = bc[2]
    @test isapprox(integral(m3, node_coords, 1.0) / expected_area, expected_val; atol=1e-10)


    # @code_warntype face_integral(node_coords,1,0,bc)
    b = @b integral($m3, $node_coords, 1.0)
    display(b)
end


