using StaticArrays
using Test
using Ju3VEM
using Ju3VEM.VEMUtils: Mesh, StandardEl, create_volume_bmat, h1_projectors!, create_node_mapping
using Chairmarks
using JET

function _create_pyramid_topology()
    topo = Topology{3}()
    base = [SA[0.0, 0.0, 0.0],
            SA[1.0, 0.0, 0.0],
            SA[1.0, 1.0, 0.0],
            SA[0.0, 1.0, 0.0]]
    apex = SA[0.5, 0.5, 1.0]

    node_ids = Int[]
    append!(node_ids, add_node!.(base, Ref(topo)))
    push!(node_ids, add_node!(apex, topo))

    base_ids = node_ids[1:4]
    apex_id = node_ids[5]

    faces = [
        base_ids,
        [base_ids[1], base_ids[2], apex_id],
        [base_ids[2], base_ids[3], apex_id],
        [base_ids[3], base_ids[4], apex_id],
        [base_ids[4], base_ids[1], apex_id]
    ]
    add_volume!(faces, topo)
    return topo
end


function _create_facedata_col(mesh::Mesh{3,StandardEl{K}}) where K
    K_MAX = max(2*K-1, 2)
    base = get_base(BaseInfo{2, K_MAX, 1}()).base
    topo = mesh.topo
    # facedata_col = Dict{Int,FaceData{3, length(base)}}()
    # for face in RootIterator{3}(topo)
    #     dΩ = precompute_face_monomials(face.id, topo, Val(K_MAX))
    #     facedata_col[face.id] = h1_projectors!(face.id, mesh, dΩ)
    # end
    facedata_col = Dict(
        face.id => h1_projectors!(
            face.id, mesh, 
            precompute_face_monomials(face.id, topo, Val(K_MAX))
            ) for face in RootIterator{3}(topo)
    )
    return facedata_col
end


using Ju3VEM.VEMUtils: create_volume_dmat
@testset "3D Mesh pipeline: creation, integration, and B-matrix" begin
    topo = _create_pyramid_topology()
    mesh = Mesh(topo, StandardEl{1}())
    base3d = get_base(BaseInfo{3,1,1}()).base

    # Mesh creation sanity checks
    @test length(get_vertices(mesh)) == 5
    # @test length(get_gauss_nodes(mesh)) == 0  # K=1 → no interior edge nodes

    # Precompute face data and volume integrals
    facedata_col = _create_facedata_col(mesh)
    vol_data = precompute_volume_monomials(1, topo, facedata_col, Val(1))

    # 3D integration check: volume should be 1/3 for unit square base and height 1 pyramid
    @test isapprox(
        compute_volume_integral_unshifted(Monomial(1.0, SA[0,0,0]), vol_data),
        1/3; atol=1e-10
    )

    # Build 3D B-matrix and validate basic properties for K=1
    ntl = create_node_mapping(1, mesh, facedata_col)
    abs_volume = vol_data.integrals[1]
    bmat3d = create_volume_bmat(1, mesh, vol_data.vol_bc, vol_data.hvol, abs_volume, facedata_col, ntl)

    dmat3d = create_volume_dmat(1, mesh, vol_data.vol_bc, vol_data.hvol, facedata_col, vol_data, abs_volume, ntl)

    gmat3d = SMatrix{length(base3d),length(base3d)}(bmat3d*dmat3d)

    proj_s = inv(gmat3d) * bmat3d


    const_uv = ones(length(mesh.nodes))
    proj_uv = dmat3d * proj_s * const_uv

    @test const_uv ≈ proj_uv
    display(proj_uv)




    lin_uv_fun(x) = x[1] + 3*x[2] + 2*x[3]

    lin_uv        = zeros(length(mesh.nodes))

    for (local_pos,node_id) in ntl.map
        lin_uv[local_pos] = lin_uv_fun(mesh[node_id])
    end

    proj_uv = dmat3d * proj_s * lin_uv

    @test lin_uv ≈ proj_uv
    display(proj_uv)
end




rep =  begin
    topo = _create_pyramid_topology()
    mesh = Mesh(topo, StandardEl{2}())
    base3d = get_base(BaseInfo{3,2,1}()).base

    # Mesh creation sanity checks

    # Precompute face data and volume integrals
    facedata_col = _create_facedata_col(mesh)
    vol_data = precompute_volume_monomials(1, topo, facedata_col, Val(2))



    # Build 3D B-matrix and validate basic properties for K=1
    ntl = create_node_mapping(1, mesh, facedata_col)
    abs_volume = vol_data.integrals[1]
    bmat3d = create_volume_bmat(1, mesh, vol_data.vol_bc, vol_data.hvol, abs_volume, facedata_col, ntl)

    dmat3d = create_volume_dmat(1, mesh, vol_data.vol_bc, vol_data.hvol, facedata_col, vol_data, abs_volume, ntl)

    gmat3d = SMatrix{length(base3d),length(base3d)}(bmat3d*dmat3d)



    proj_s = inv(gmat3d) * bmat3d

    proj = dmat3d * proj_s
    proj_s2, proj2 = create_volume_vem_projectors(1, mesh, vol_data, facedata_col, ntl)

    @test proj ≈ proj2
    @test proj_s ≈ proj_s2

    uv = rand(length(mesh.nodes))

    proj_uv = proj * uv

    proj2_uv = proj * proj_uv

    @test proj_uv ≈ proj2_uv

    uv2 = dmat3d[:,end] + dmat3d[:,end-1]

    proj_uv2 = proj2 * uv2
    @test uv2 ≈ proj_uv2

    uv3 = dmat3d[:,1]
    
    proj_uv3 = proj * uv3
    @test uv3 ≈ proj_uv3

    b = @b create_volume_bmat(1, $mesh, $vol_data.vol_bc, $vol_data.hvol, $abs_volume, $facedata_col, $ntl)
    display(b)

    b2 = @b create_volume_dmat(1, $mesh, $vol_data.vol_bc, $vol_data.hvol, $facedata_col, $vol_data, $abs_volume, $ntl)
    display(b2)


end


