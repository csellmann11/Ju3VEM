# Run with an environment that develops Ju3VEM, for example:
# julia --project=../ToOpt3 --startup-file=no tests/projector_regression.jl
using Test, Ju3VEM, LinearAlgebra, StaticArrays
import Ju3VEM.VEMUtils as VU

function test_projector_reproduction(mesh)
    cv = CellValues{3}(mesh)
    for fd in values(cv.facedata_col)
        D = Matrix(VU.create_D_mat(mesh,fd))
        @test Matrix(fd.ΠsL2)*D ≈ I atol=2e-10
    end
    for element in RootIterator{4}(mesh.topo)
        reinit!(element.id,cv)
        D = Matrix(VU.create_volume_dmat(element.id,mesh,
            cv.facedata_col,cv.volume_data,cv.vnm))
        Ps,P = create_volume_vem_projectors(element.id,mesh,
            cv.volume_data,cv.facedata_col,cv.vnm)
        @test Matrix(Ps)*D ≈ I atol=2e-10
        @test Matrix(P)*D ≈ D atol=2e-10
        @test Matrix(P)^2 ≈ Matrix(P) atol=2e-10
    end
end

@testset "Production VEM projectors reproduce polynomials" begin
    # A unit cube used to give diag(1,1/9,1/9,1/9) instead of identity
    # through the Octavian.matmul! path with an immutable static inverse.
    @testset "Unit cube" begin
        test_projector_reproduction(create_rectangular_mesh(1,1,1,1.,1.,1.,StandardEl{1}))
    end
    @testset "Scaled anisotropic cells" begin
        test_projector_reproduction(create_rectangular_mesh(2,2,2,0.4,0.2,0.1,StandardEl{1}))
    end
    @testset "Hanging nodes" begin
        mesh = create_rectangular_mesh(2,2,2,1.,1.,1.,StandardEl{1})
        Ju3VEM.VEMGeo._refine!(first(RootIterator{4}(mesh.topo)),mesh.topo)
        test_projector_reproduction(Mesh(mesh.topo,StandardEl{1}()))
    end
    @testset "Quadratic cube" begin
        test_projector_reproduction(create_rectangular_mesh(1,1,1,1.,1.,1.,StandardEl{2}))
    end
end
