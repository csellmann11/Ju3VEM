using Ju3VEM
using Ju3VEM.VEMGeo: get_base_len,get_plane_parameters,face2d_sym_integral
using StaticArrays
using LinearAlgebra
using Symbolics, Chairmarks
using Test
using FixedSizeArrays, Bumper
using Ju3VEM.VEMGeo: D2FaceParametrization, project_to_2d_abs, project_to_2d_rel, project_to_3d, project_to_3d_flat

@testset "2D Integration Tests" begin
    node_coords = 2 .* [SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0), SVector(0.0, 1.0)]

    m1 = Monomial(1.0, SVector(0, 0))

    expected_area = 4.0
    fd = Ju3VEM.VEMGeo.precompute_face_monomials(node_coords, Val(1))
    @test isapprox(get_area(fd), expected_area; atol=1e-10)


    m2 = Monomial(1.0, SVector(1, 0))
    bc = SVector(1.0, 1.0)
    # linear monomial integrates to 0 when centered at face barycenter
    @test isapprox(Ju3VEM.VEMGeo.compute_face_integral(m2, fd, bc), 0.0; atol=1e-10)


    expected_val = bc[1]
    m3 = Monomial(1.0, SVector(1, 0))
    @test isapprox(Ju3VEM.VEMGeo.compute_face_integral(m3, fd, SA[0.0,0.0],1.0) / expected_area, expected_val; atol=1e-10)
end


@testset "Test Plane Parameterization" begin
    quad_nodes = [SA[0.0,0.0,0.0], SA[1.0,0.0,0.0], SA[1.0,1.0,1.0], SA[0.0,1.0,1.0]] .* 2 .+ Scalar(SA[1.0,1.0,1.0])
    hf = max_node_distance(quad_nodes)
    
    plane = D2FaceParametrization(quad_nodes)

    quad_nodes2d = project_to_2d_abs.(quad_nodes,Ref(plane))
    quad_nodes3d = project_to_3d.(quad_nodes2d,Ref(plane))
    @test quad_nodes3d ≈ quad_nodes
    
    quad_nodes2d_rel = project_to_2d_rel.(quad_nodes,Ref(plane))
    quad_nodes3d_rel = project_to_3d.(quad_nodes2d_rel,Ref(plane))

    @test quad_nodes3d_rel ≈ quad_nodes .- Scalar(plane.p0)
end



@testset "Test tilted face integration" begin

    quad_nodes = [SA[0.0,0.0,0.0], SA[1.0,0.0,0.0], SA[1.0,1.0,1.0], SA[0.0,1.0,1.0]] .* 2 .+ Scalar(SA[1.0,1.0,1.0])
    hf = max_node_distance(quad_nodes)

    plane = D2FaceParametrization(quad_nodes)
    a   = face2d_sym_integral(quad_nodes,0,0,SA[0.0,0.0],plane)
    bcu = face2d_sym_integral(quad_nodes,1,0,SA[0.0,0.0],plane)/a
    bcv = face2d_sym_integral(quad_nodes,0,1,SA[0.0,0.0],plane)/a

    face_data = precompute_face_monomials(quad_nodes,Val(5))



    # test if integration over the projected faces gives the same results
    quad_nodes_2d = project_to_2d_abs.(quad_nodes,Ref(plane))

    plane2d = D2FaceParametrization(quad_nodes_2d)

    a2d   = face2d_sym_integral(quad_nodes_2d,0,0,SA[0.0,0.0],plane2d)
    bcu2d = face2d_sym_integral(quad_nodes_2d,1,0,SA[0.0,0.0],plane2d)/a2d
    bcv2d = face2d_sym_integral(quad_nodes_2d,0,1,SA[0.0,0.0],plane2d)/a2d

    @test bcu ≈ bcu2d
    @test bcv ≈ bcv2d  

    bc3d = project_to_3d(face_data.bc,plane)
    @test bc3d ≈ SA[2.0,2.0,2.0]



    true_area = (quad_nodes_2d[2][1] - quad_nodes_2d[1][1])*(quad_nodes_2d[3][2] - quad_nodes_2d[1][2])

    @test true_area ≈ a


    # test monomial over 3d face
    m3d = Monomial(1.0,SA[2,1,1])
    Int = compute_face_integral(m3d,face_data)



    ref = 4*sqrt(2)/(9*hf^4)

    @test Int ≈ ref

    # test allocations
    allocs = @allocated compute_face_integral(m3d,face_data)
    @test allocs == 0

end




@testset "Test horizontal elevated face integration" begin

    quad_nodes2 = [SA[0.0,0.0,0.0], SA[1.0,0.0,0.0], SA[1.0,1.0,0.0], SA[0.0,1.0,0.0]] .* 2 .+ Scalar(SA[1.0,1.0,1.0])
    hf = max_node_distance(quad_nodes2)
    u, v, n, nus,p0 = Ju3VEM.VEMGeo.get_plane_parameters(quad_nodes2)

    face_data2 = precompute_face_monomials(quad_nodes2,Val(5))
    bc3d = project_to_3d(face_data2.bc,face_data2.plane)
    bc3d_plane  = project_to_3d_flat(face_data2.bc,face_data2.plane)


    m3d = Monomial(1.0,SA[0,1,1])

    Int = compute_face_integral(m3d,face_data2,SA[0.0,0.0,0.0],2.0)

    @test Int ≈ 2.0

    quad_nodes3 = ([SA[0.0,0.0,0.0], SA[1.0,0.0,0.0], SA[1.0,1.0,0.0], SA[0.0,1.0,0.0]] .* 2 .+ 2*Scalar(SA[1.0,1.0,0.0])) 
    face_data3 = precompute_face_monomials(quad_nodes3,Val(1))


    m3d = Monomial(1.0,SA[0,1,0])
    m2d = Monomial(1.0,SA[0,1])

    Int  = compute_face_integral(m3d,face_data3,SA[0.0,0.0,0.0],2.0)
    I2 = compute_face_integral(m2d,face_data3,SA[0.0,0.0],2.0)


    allocs = @allocated compute_face_integral(m2d,face_data3,SA[0.0,0.0],2.0)
    @test allocs == 0


    @test Int  ≈ 6.0 
    @test I2 ≈ 6.0 


    I3 = compute_face_integral(m2d,face_data3,face_data3.bc,2.0)
    I4 = compute_face_integral_unshifted(m2d,face_data3)


    allocs = @allocated compute_face_integral_unshifted(m2d,face_data3)
    @test allocs == 0
end







