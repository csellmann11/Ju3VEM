using StaticArrays
using LinearAlgebra
using Test
using Chairmarks
using LoopVectorization
using Bumper, OrderedCollections, SmallCollections
using Ju3VEM
# import .Triangulation3D: FaceTriangulations3D

# @testset "3D Face Integration (planar polygon)" begin
#     # Build a rectangle in a random oriented plane
#     o = SA[0.3, -0.2, 0.4]
#     # Orthonormal basis
#     n = normalize(SA[0.2, 1.0, -0.5])
#     tmp = abs(n[1]) < 0.9 ? SA[1.0,0.0,0.0] : SA[0.0,1.0,0.0]
#     u = normalize(cross(n, tmp))
#     v = cross(n, u)
#     a = 1.2
#     b = 0.7
#     rect2 = [SA[-a,-b], SA[a,-b], SA[a,b], SA[-a,b]]
#     rect3 = [o + p[1]*u + p[2]*v for p in rect2]

#     # Constant monomial integrates to area
#     m0 = Monomial(1.0, SVector(0,0,0))
#     I0 = integral(m0, rect3)
#     A = Triangulation3D.polygon_area3D(rect3)
#     @test isapprox(I0, A; rtol=1e-12, atol=1e-12)

#     # Linear monomial centered at rectangle center integrates to zero
#     bc = o
#     m1x = Monomial(1.0, SVector(1,0,0))
#     m1y = Monomial(1.0, SVector(0,1,0))
#     m1z = Monomial(1.0, SVector(0,0,1))
#     @test isapprox(integral(m1x, rect3, 1.0, bc), 0.0; atol=1e-12)
#     @test isapprox(integral(m1y, rect3, 1.0, bc), 0.0; atol=1e-12)
#     @test isapprox(integral(m1z, rect3, 1.0, bc), 0.0; atol=1e-12)

#     # quick perf sanity
#     _ = integral(m1x, rect3, 1.0, bc)
# end


reps = []
@testset "3D Volume Integration (polyhedron via face triangles)" begin
    # Unit cube
    pts2 = [
        SA[0.0,0.0,0.0], SA[1.0,0.0,0.0], SA[1.0,1.0,0.0], SA[0.0,1.0,0.0],
        SA[0.0,0.0,1.0], SA[1.0,0.0,1.0], SA[1.0,1.0,1.0], SA[0.0,1.0,1.0]
    ] .* 2.0
    faces = [
        [1,2,3,4], # bottom z=0
        [5,6,7,8], # top z=1
        [1,2,6,5],
        [2,3,7,6],
        [3,4,8,7],
        [4,1,5,8],
    ]

    # Build topology and precompute face triangulations
    topo2 = Topology{3}()
    add_node!.(pts2, Ref(topo2))
    add_volume!(faces, topo2)
    vol_id2 = 1
    ft2 = FaceTriangulations3D(topo2)

    # Constant: volume = 1
    m0 = Monomial(1.0, SVector(0,0,0))
    V0 = integrate_polynomial_over_volume(m0, vol_id2, topo2, ft2)
    expected_vol = 8.0
    @test isapprox(V0, 8.0; rtol=1e-12, atol=1e-12)

    # Linear: ∫ x dV = 1/2, ∫ y dV = 1/2, ∫ z dV = 1/2
    mx = Monomial(1.0, SVector(1,0,0))
    my = Monomial(1.0, SVector(0,1,0))
    mz = Monomial(1.0, SVector(0,0,1))
    @test isapprox(integrate_polynomial_over_volume(mx, vol_id2, topo2, ft2)/expected_vol, 1.0; rtol=1e-12, atol=1e-12)
    @test isapprox(integrate_polynomial_over_volume(my, vol_id2, topo2, ft2)/expected_vol, 1.0; rtol=1e-12, atol=1e-12)
    @test isapprox(integrate_polynomial_over_volume(mz, vol_id2, topo2, ft2)/expected_vol, 1.0; rtol=1e-12, atol=1e-12)

   
    pts = [
        SA[0.0,0.0,0.0], SA[1.0,0.0,0.0], SA[1.0,1.0,0.0], SA[0.0,1.0,0.0],
        SA[0.0,0.0,1.0], SA[1.0,0.0,1.0], SA[1.0,1.0,1.0], SA[0.0,1.0,1.0]
    ] 

    topo = Topology{3}()
    add_node!.(pts, Ref(topo))
    add_volume!(faces, topo)
    vol_id = 1
    ft = FaceTriangulations3D(topo)

    # Quadratic: ∫ x^2 dV = 1/3, cross term ∫ x*y dV = 1/4
    mxx = Monomial(1.0, SVector(2,0,0))
    mxy = Monomial(1.0, SVector(1,1,0))
    @test isapprox(integrate_polynomial_over_volume(mxx, vol_id, topo, ft), 1/3; rtol=1e-12, atol=1e-12)
    @test isapprox(integrate_polynomial_over_volume(mxy, vol_id, topo, ft), 1/4; rtol=1e-12, atol=1e-12)

end

