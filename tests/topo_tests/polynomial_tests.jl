using StaticArrays
using Chairmarks
using Test
using Symbolics, LinearAlgebra, Statistics
using Ju3VEM
using FixedSizeArrays
using Ju3VEM.VEMGeo: poly_pow

@testset "Polynomial Basic Operations" begin
    @testset "Polynomial Construction and Evaluation" begin
        bi = BaseInfo{2,2,1}()
        base = get_base(bi)
        p = Polynomial(ones(length(base)), base.base)

        x = SA[1.0, 2.0]
        expected_val = 1 + x[1] + x[2] + x[1]^2 + x[2]^2 + x[1]*x[2]
        @test p(x) ≈ expected_val
    end

    @testset "Polynomial Gradients" begin
        bi = BaseInfo{2,2,1}()
        base = get_base(bi)
        p = Polynomial(ones(length(base)), base.base)

        x = SA[1.0, 2.0]
        expected_grad = [SA[1.0, 2.0, 1.0], SA[1.0, 1.0, 2.0]]
        @test ∇p(p, 1.0)[1].coeffs ≈ expected_grad[1]
        @test ∇p(p, 1.0)[2].coeffs ≈ expected_grad[2]

        expected_eval = SA[1.0 + 2.0*x[1] + 1.0*x[2], 1.0 + 1.0*x[1] + 2.0*x[2]]
        @test ∇p(p, 1.0)(x) ≈ expected_eval
    end

    @testset "Polynomial Matrix Operations" begin
        bi = BaseInfo{2,2,1}()
        base = get_base(bi)
        x = SA[1.0, 2.0]
        p = Polynomial(ones(length(base)), base.base)

        pmat = SA[p p; p p]
        @test ∇p(pmat, 1.0)[1,1,1] == ∇p(p, 1.0)[1]
        @test ∇x(p, 1.0, x) ≈ ∇p(p, 1.0)(x)
    end
end

@testset "Polynomial Arithmetic" begin
    @testset "Polynomial Multiplication" begin
        bi = BaseInfo{2,1,1}()
        base = get_base(bi)
        coeffs = (1:length(base)) .|> Float64
        p = Polynomial(coeffs, base.base)

        x = SA[1.0, 2.0]
        true_coeffs = [1.0, 4.0, 6.0, 4.0, 12.0, 9.0]
        @test (p*p).coeffs ≈ true_coeffs
        @test p(x)*p(x) ≈ (p*p)(x)
    end

    @testset "Polynomial Power Operations" begin
        bi = BaseInfo{2,1,1}()
        base = get_base(bi)
        coeffs = (1:length(base)) .|> Float64
        p = Polynomial(coeffs, base.base)

        @test poly_pow(p, Val(3)) ≈ p*p*p
        
        # Test additional power cases
        @test poly_pow(p, Val(1)) ≈ p
        @test poly_pow(p, Val(2)) ≈ p*p
    end
end

@testset "Coordinate Transformations" begin
    # @testset "2D to 3D Transformation" begin
    #     quad_nodes = [SA[0.0,0.0,0.0], SA[1.0,0.0,0.0], SA[1.0,1.0,1.0], SA[0.0,1.0,1.0]]
    #     u, v, n, nus, p0 = Ju3VEM.VEMGeo.get_plane_parameters(quad_nodes)

    #     # Test 2d to 3d
    #     m2d = Monomial(1.0, SA[1,1])
    #     coeffs = compute_transformation_coeffs2d_to_3d(m2d, p0, u, v)
        
    #     # Benchmark the transformation
    #     be = @be compute_transformation_coeffs2d_to_3d($m2d, $p0, $u, $v)
    #     display(be)

    #     lin3d_base = get_base(BaseInfo{3,1,1}()).base
    #     p1 = Polynomial(vcat(0.0, u), lin3d_base)
    #     p2 = Polynomial(vcat(0.0, v), lin3d_base)
    #     p = p1 * p2

    #     @test p.coeffs ≈ coeffs

    #     x2d = SA[0.5, 0.3]
    #     val2d = m2d(x2d)
    #     x3d = x2d[1]*u + x2d[2]*v
    #     val3d = p(x3d)

    #     @test val2d ≈ val3d
    # end

    @testset "3D to 2D Transformation" begin
        quad_nodes = [SA[0.0,0.0,0.0], SA[1.0,0.0,0.0], SA[1.0,1.0,1.0], SA[0.0,1.0,1.0]]
        u, v, n, nus, p0 = Ju3VEM.VEMGeo.get_plane_parameters(quad_nodes)

        # Test 3d to 2d
        m3d = Monomial(1.0, SA[1,1,1])

        lin2d_base = get_base(BaseInfo{2,1,1}()).base
        p1 = Polynomial(SA[p0[1], u[1], v[1]], lin2d_base)
        p2 = Polynomial(SA[p0[2], u[2], v[2]], lin2d_base)
        p3 = Polynomial(SA[p0[3], u[3], v[3]], lin2d_base)

        p = p1 * p2 * p3
        coeffs = compute_transformation_coeffs3d_to_2d(m3d, p0, u, v)

        b = @b compute_transformation_coeffs3d_to_2d!($coeffs,$m3d, $p0, $u, $v)
        display(b)


        # Benchmark the transformation
        be = @b compute_transformation_coeffs3d_to_2d($m3d, $p0, $u, $v)
        display(be)

        @test p.coeffs ≈ coeffs

        base_cub_2d = get_base(BaseInfo{2,3,1}()).base
        p2d = Polynomial(coeffs, base_cub_2d)

        x2d = SA[0.5, 0.3]
        x3d = x2d[1]*u + x2d[2]*v
        val3d = m3d(x3d)
        val2d = p2d(x2d)
        
        @test val2d ≈ val3d
    end
end

@testset "Edge Cases and Error Handling" begin
    @testset "Zero Polynomial" begin
        bi = BaseInfo{2,1,1}()
        base = get_base(bi)
        zero_poly = Polynomial(zeros(length(base)), base.base)
        
        x = SA[0.5, 0.3]
        @test zero_poly(x) ≈ 0.0 
        @test (zero_poly * zero_poly)(x) ≈ 0.0
    end

    @testset "Constant Polynomial" begin
        bi = BaseInfo{2,1,1}()
        base = get_base(bi)
        const_poly = Polynomial([5.0, 0.0, 0.0], base.base)
        
        x = SA[0.5, 0.3]
        @test const_poly(x) ≈ 5.0
        @test ∇p(const_poly, 1.0)(x) ≈ SA[0.0, 0.0]
    end

    @testset "Polynomial Evaluation at Origin" begin
        bi = BaseInfo{2,2,1}()
        base = get_base(bi)
        p = Polynomial(ones(length(base)), base.base)
        
        origin = SA[0.0, 0.0]
        @test p(origin) ≈ 1.0  # Only constant term contributes
    end
end




