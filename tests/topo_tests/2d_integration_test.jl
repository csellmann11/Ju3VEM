using StaticArrays, WriteVTK 
using OrderedCollections, Bumper
using LinearAlgebra, Statistics
using SmallCollections, Chairmarks
using LoopVectorization
using Test
using Ju3VEM

@testset "2D Integration Tests" begin
    node_coords = 2 .* [SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0), SVector(0.0, 1.0)]

    m1 = Monomial(1.0, SVector(0, 0))

    expected_area = 4.0
    fd = Ju3VEM.VEMGeo.precompute_face_monomials(node_coords, Val(1))
    @test isapprox(get_area(fd), expected_area; atol=1e-10)


    m2 = Monomial(1.0, SVector(1, 0))
    bc = SVector(1.0, 1.0)
    # linear monomial integrates to 0 when centered at face barycenter
    @test isapprox(Ju3VEM.VEMGeo.compute_face_integral(m2, fd), 0.0; atol=1e-10)

    expected_val = 0.0
    @test isapprox(Ju3VEM.VEMGeo.compute_face_integral(m2, fd, SA[0.0,0.0]) / expected_area, expected_val; atol=1e-10)


    expected_val = bc[1]
    @test isapprox(Ju3VEM.VEMGeo.compute_face_integral(m3, fd, SA[0.0,0.0]) / expected_area, expected_val; atol=1e-10)

    m3 = Monomial(1.0, SVector(0, 1))
    expected_val = bc[2]
    @test isapprox(expected_val, expected_val; atol=1e-10)


    # @code_warntype face_integral(node_coords,1,0,bc)
    # removed benchmark of integral
end


