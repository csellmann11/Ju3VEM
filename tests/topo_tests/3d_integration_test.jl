using StaticArrays
using LinearAlgebra
using Test
using Chairmarks
using LoopVectorization
using Bumper, OrderedCollections, SmallCollections
using Ju3VEM
using Ju3VEM.VEMUtils: _create_facedata_col
# import .Triangulation3D: FaceTriangulations3D


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

    bc = SA[0.5,0.5,0.5]
    # topo_shifted = deepcopy(topo)
    # for i in eachindex(topo_shifted.nodes)
    #     topo_shifted.nodes[i] = Node(topo_shifted.nodes[i] - bc)
    # end


    ft = FaceTriangulations3D(topo)

    # Quadratic: ∫ x^2 dV = 1/3, cross term ∫ x*y dV = 1/4
    mxx = Monomial(1.0, SVector(2,0,0))
    mxy = Monomial(1.0, SVector(1,1,0))
    @test isapprox(integrate_polynomial_over_volume(mxx, vol_id, topo, ft), 1/3; rtol=1e-12, atol=1e-12)
    @test isapprox(integrate_polynomial_over_volume(mxy, vol_id, topo, ft), 1/4; rtol=1e-12, atol=1e-12)





    mesh = Mesh(topo, StandardEl{2}())
    # cv = CellValues{1}(mesh)
    fdc = _create_facedata_col(mesh)
    vol_data = precompute_volume_monomials(1, mesh.topo, fdc, Val(2), false)
    # @show vol_data.vol_bc
    
    mi = Monomial(1.0, SVector(1,1,0))
    isym = compute_volume_integral_unshifted(mi, vol_data,1.0)
    inum = integrate_polynomial_over_volume(mi, vol_id, topo, ft)
    @show isym
    @show inum
    @show inum/isym
    @show isym/inum
    # @test isapprox(isym, inum; rtol=1e-12, atol=1e-12)

end

