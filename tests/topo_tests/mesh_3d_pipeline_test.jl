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

    dmat3d = create_volume_dmat(1, mesh,  facedata_col, vol_data, ntl)

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


using Ju3VEM.VEMUtils: get_local_id, static_matmul, create_D_mat, create_B_mat

# INFO: Testing second order elements
rep =  begin
    topo = _create_pyramid_topology()
    mesh = Mesh(topo, StandardEl{2}())
    # mesh = create_unit_rectangular_mesh(1,1,1, StandardEl{2}) 
    # topo = mesh.topo
    base3d = get_base(BaseInfo{3,2,1}()).base

    base2d = get_base(BaseInfo{2,2,1}()).base

    facedata_col = _create_facedata_col(mesh)


    vol_data = precompute_volume_monomials(1, mesh.topo, facedata_col, Val(2))

    bc_vol = vol_data.vol_bc
    hvol = vol_data.hvol


    # Build 3D B-matrix and validate basic properties for K=1
    ntl = create_node_mapping(1, mesh, facedata_col)
    abs_volume = vol_data.integrals[1]
    # bmat3d = create_volume_bmat(1, mesh, vol_data.vol_bc, vol_data.hvol, abs_volume, facedata_col, ntl)

    # dmat3d = create_volume_dmat(1, mesh, facedata_col, vol_data, ntl)

    # gmat3d = SMatrix{length(base3d),length(base3d)}(bmat3d*dmat3d)


    proj_s, proj = create_volume_vem_projectors(1, mesh, vol_data, facedata_col, ntl)


    @show sum(proj,dims = 2)

    # INFO: test of 2d face projectors 
end



# Info: Second order for linear recovery
begin 

    u_true = [
    1.0
    2.0
    2.0
    3.0000000000000004
    2.0
    3.0000000000000004
    3.0000000000000004
    4.0
    1.5000000000000002
    2.5
    2.5
    1.5000000000000002
    2.5
    3.5000000000000004
    3.5000000000000004
    2.5
    2.5
    1.5000000000000002
    3.5000000000000004
    2.5
    2.0
    3.0000000000000004
    2.0
    3.0000000000000004
    2.0
    3.0000000000000004
    2.5 ]


    mesh = create_unit_rectangular_mesh(1,1,1, StandardEl{2}) 

    # Precompute face data and volume integrals
    facedata_col = _create_facedata_col(mesh)



    vol_data = precompute_volume_monomials(1, mesh.topo, facedata_col, Val(2))
    vnm = create_node_mapping(1, mesh, facedata_col)
    u_local_perm = zero(u_true)

    local_mom_ids = Int[]
    for (node_id,local_pos) in vnm.map
        u_local_perm[local_pos] = u_true[node_id]
        if node_id > 20 
            push!(local_mom_ids,local_pos)
        end
    end

    hvol = vol_data.hvol
    bcvol = vol_data.vol_bc
    abs_volume = vol_data.integrals[1]

    proj_s, proj = create_volume_vem_projectors(1, mesh, vol_data, facedata_col, vnm)

    coeffs = proj_s *  u_local_perm 
    display(coeffs)

    u_proj = proj * u_local_perm
    display(u_local_perm)
    display(u_proj)

    dm3d = create_volume_dmat(1, mesh, facedata_col, vol_data, vnm)


    # test_fun(x) = 1/hvol*(x[1]-bcvol[1]) 
    test_fun(x) = 1 + x[1] + x[2] + x[3]
    coeffs = zeros(10)
    coeffs[1] = sum(bcvol) + 1.0
    coeffs[2] = hvol 
    coeffs[3] = hvol 
    coeffs[4] = hvol 


    uvals = zeros(length(mesh.nodes))
    for (node_id,local_pos) in vnm.map
        uvals[local_pos] = test_fun(mesh[node_id])

        if node_id > 20 
            uvals[local_pos] = dm3d[local_pos,:] ⋅ coeffs
        end

        test_val = dm3d[local_pos,:] ⋅ coeffs
        @test test_val ≈ uvals[local_pos] atol=1e-10
    end

    computed_coeffs = proj_s * uvals
    uvals_proj = proj * uvals
    @test uvals ≈ uvals_proj atol=1e-10
    @test computed_coeffs ≈ coeffs atol=1e-10
    @test u_local_perm[local_mom_ids] ≈ uvals_proj[local_mom_ids] atol=1e-10
end




# Info: Second order for quad recovery

begin 
    mesh = create_unit_rectangular_mesh(1,1,1, StandardEl{2}) 

    base3d = get_base(BaseInfo{3,2,1}()).base
    # Precompute face data and volume integrals
    facedata_col = _create_facedata_col(mesh)



    vol_data = precompute_volume_monomials(1, mesh.topo, facedata_col, Val(2))
    vnm = create_node_mapping(1, mesh, facedata_col)

    local_mom_ids = Int[]
    for (node_id,local_pos) in vnm.map
        if node_id > 20 
            push!(local_mom_ids,local_pos)
        end
    end

    hvol = vol_data.hvol
    bcvol = vol_data.vol_bc
    abs_volume = vol_data.integrals[1]

    proj_s, proj = create_volume_vem_projectors(1, mesh, vol_data, facedata_col, vnm)

    dm3d = create_volume_dmat(1, mesh, facedata_col, vol_data, vnm)

    ft = FaceTriangulations3D(mesh.topo)

    # test_fun(x) = 1/hvol*(x[1]-bcvol[1]) 
    test_fun(x) = x[1]^2
    poly_unscaled_coeffs = zeros(10)
    poly_unscaled_coeffs[5] = 1.0
    poly_unscaled = Polynomial(poly_unscaled_coeffs, base3d)


    coeffs = zeros(10)
    coeffs[1] = bcvol[1]^2
    coeffs[2] = 2*bcvol[1]*hvol
    coeffs[5] = hvol^2

    poly = Polynomial(coeffs, base3d)

    u_sol = zeros(length(mesh.nodes))
    for (node_id,local_pos) in vnm.map
        node = mesh[node_id]
        u_sol[local_pos] = test_fun(node)

        scaled_loc = (node - bcvol)/hvol
        @test u_sol[local_pos] ≈ poly(scaled_loc) atol=1e-10
    end

    for (face_id,face_data) in facedata_col 

        face_moment_id = face_data.face_node_ids[end]
        local_moment_id = vnm.map[face_moment_id]
        u_sol[local_moment_id] = 
            integrate_polynomial_over_face(poly_unscaled,face_id,mesh.topo,ft)[1]
    end
    u_sol[end] = integrate_polynomial_over_volume(poly_unscaled,1,mesh.topo,ft)

    u_proj = proj * u_sol 
    @test u_sol ≈ u_proj atol=1e-10

    proj_coeffs = proj_s * u_sol
    @test proj_coeffs ≈ coeffs atol=1e-10

    display(u_sol)

    u_global_perm = zero(u_sol)
    for (node_id,local_pos) in vnm.map
        u_global_perm[node_id] = u_sol[local_pos]
    end

    display(u_global_perm)

end