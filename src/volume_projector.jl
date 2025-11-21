using ..VEMGeo: project_to_3d_flat, compute_face_integral_coeffs!



function precompute_base_projection_coeffs(::BaseInfo{3,O},
    dΩ::FaceIntegralData,
    bc_vol::SVector{3,Float64},
    hvol::Float64) where {O}


    base    = get_base(BaseInfo{3,O,1}()).base 
    base2d  = get_base(BaseInfo{2,O,1}()).base 

    coeff_list = 
            MVector{length(base),SVector{length(base2d),Float64}}(undef)

    @no_escape begin 
        vec = @alloc(Float64,length(base2d))
        for (i,m) in enumerate(base)
            vec .= 0.0 
            compute_face_integral_coeffs!(vec,m,dΩ,bc_vol)

            order = sum(m.exp)

            coeff_list[i] = SVector{length(base2d),Float64}(vec)/hvol^order
        end
    end

    return SVector(coeff_list)
end

function get_dot_product_projection_coeffs(∇m::SVector{3,Monomial{Float64,3}}, 
    n::SVector{3,Float64},
    coeff_list::SVector{L,<:SVector}) where {L} 


    coeffs = zero(eltype(coeff_list))

    for (i,∇mi) in enumerate(∇m)
        idx = get_exp_to_idx_dict(∇mi.exp)
        coeffs += coeff_list[idx]*n[i]*∇mi.val
    end 

    return coeffs
end


function compute_face_integral_m2d_time_3dcoeffs(m2d::Monomial{T,2},
    fd::FaceIntegralData{D,L},coeffs::SVector,
    hf::T = fd.hf) where {D,T,L}

    # function computes ∫ m2d((x-bc_f)/hf) * coeffs_of(nx*dm3dx + ny*dm3dy + nz*dm3dz)
    # the coeffs already contain the ∇m3d_dot_n and also the h scaling and the correct offset of the 
    # volume bary center


    order2 = sum(m2d.exp)
    full_base = get_base(BaseInfo{2,6,1}()).base

    λ = hf^order2
    ∫m = zero(T)
    for i in eachindex(coeffs) 
        ci = coeffs[i]
        new_exp   = m2d.exp + full_base[i].exp  
        new_coeff = m2d.val*ci 
 
        idx = get_exp_to_idx_dict(new_exp)
        ∫m += new_coeff*fd.integrals[idx]/λ
    end
    return ∫m
end

function create_volume_bmat(volume_id::Int,
    mesh::Mesh{3,ET},
    bcvol,hvol,abs_volume,
    facedata_col::Dict{Int,<:FaceData{3,L}},
    volume_node_mapping::NodeID2LocalID) where {K,L,ET<:ElType{K}}


    topo = mesh.topo

    base_2d = get_base(BaseInfo{2,K,1}()).base
    base_3d = get_base(BaseInfo{3,K,1}()).base

 
    n_nodes = get_n_nodes(volume_node_mapping)
    moment_idx_start_minus1 = sum(volume_node_mapping.counter[i] for i in 1:3)


    bmat = FixedSizeMatrix{Float64}(undef,length(base_3d),n_nodes)
    bmat .= 0.0 

    if K == 1 
        bmat[1,:] .= 1/n_nodes
    else 
        bmat[1,moment_idx_start_minus1+1] = 1.0
    end

    iterate_volume_areas(facedata_col,topo,volume_id) do _, face_data, _
   
        face_normal = get_outward_normal(bcvol,face_data)
        ΠsL2                    = face_data.ΠsL2
        face_node_ids           = face_data.face_node_ids

        hf = face_data.dΩ.hf

        coeff_list = precompute_base_projection_coeffs(BaseInfo{3,K-1,1}(),
                        face_data.dΩ,bcvol,hvol)

        integrated_vals = MVector{length(base_2d),Float64}(undef)
        for (i,m3d) in enumerate(base_3d)
            i == 1 && continue # gradient of constant is zero
            ∇m3d = ∇(m3d,hvol) 

            # the first loop computes ∫ m2d * (∇m3d ⋅ n) ∂Ω
            # for this (∇m3d ⋅ n) is transformed into a 2d polynomial defined on the face
            for (j,m2d) in enumerate(base_2d)


                ∇m3d_dot_n = get_dot_product_projection_coeffs(
                    ∇m3d,face_normal,coeff_list)


                integrated_vals[j] = compute_face_integral_m2d_time_3dcoeffs(
                    m2d,face_data.dΩ,∇m3d_dot_n,hf)
            end
           

            # takes the res of the prev loop v = ∫ m2d * (∇m3d ⋅ n) ∂Ω for all m2d in base_2d
            # computes ΠsL2 * v and directly stores the result in bmat
            for nc in axes(ΠsL2,2)
                colv = SVector{length(base_2d)}(ΠsL2[i,nc] for i in axes(ΠsL2,1))

                col_id          = get_local_id(volume_node_mapping,face_node_ids[nc])
                bmat[i,col_id] += colv ⋅ integrated_vals
            end

            
        end
    end

    ## volume moments
    for i in 5:length(base_3d)
        m3d = base_3d[i]
        # calculates ∫ Δm * u dΩ = [f1,f2,f3] * ∫ [m,xx m,yy m,zz] * u dΩ
        for d in 1:3 
            ddm     = ∂(m3d,hvol,d,d)
            ddm.val == zero(ddm.val) && continue
            idx     = get_exp_to_idx_dict(ddm.exp)
            bmat[i,moment_idx_start_minus1+idx] -= ddm.val*abs_volume 
        end
    end
   
    return bmat
end


# function manual_K2_gmat(vol_data::VolumeIntegralData,hvol::Float64,abs_volume::Float64)

#     base3d = get_base(BaseInfo{3,2,1}()).base
    

#     gmat = MMatrix{length(base3d),length(base3d),Float64}(undef)

#     for (i,mi) in enumerate(base3d)
#         i == 1 && continue
#         ∇mi = ∇(mi,hvol)
#         for (j,mj) in enumerate(base3d)
#             j == 1 && continue
#             ∇mj = ∇(mj,hvol)

#             gmat[i,j] = sum(
#                 compute_volume_integral_unshifted(∇mi[k]*∇mj[k],vol_data,hvol) for k in 1:3
#             )
#         end
#     end

#     for (i,mi) in enumerate(base3d)
#         gmat[1,i] = compute_volume_integral_unshifted(mi,vol_data,hvol)/abs_volume
#     end
#     return SMatrix(gmat)
    

# end


function _is_vertex_like_node(id::Int,mesh::Mesh)::Bool 

    n_vertices    = get_vertices(mesh)    |> length
    n_gauss_nodes = get_gauss_nodes(mesh) |> length 

    return id <= (n_vertices + n_gauss_nodes)
end 


function _handle_face_moments!(
    dmat::AbstractMatrix,
    mesh::Mesh{3,ET},
    bc_vol::SVector{3,Float64},
    facedata_col::Dict{Int,<:FaceData},
    ntl,
    vol_id,
    hvol
) where {K,ET<:ElType{K}}


    topo = mesh.topo
    mombase2d         = get_base(BaseInfo{2,K-2,1}())
    length(mombase2d) == 0 && return nothing

    iterate_volume_areas(facedata_col,topo,vol_id) do face, facedata, _

        face_moment_ids = mesh.int_coords_connect[3][face.id]

        # gets for all 3d monomials the coeffs when expressed in 2d face basis
        face_coeffs = precompute_base_projection_coeffs(
            BaseInfo{3,K,1}(),facedata.dΩ,bc_vol,hvol
        )

        # the dofs are now \int m3di * m2dj dΓ for all m3di and m2dj
        for (i,coeffs3d) in enumerate(face_coeffs)
            for (moment_id,mi) in zip(face_moment_ids,mombase2d)

                local_pos = get_local_id(ntl,moment_id)


                # ∫ pi(moment_j(x)) * m2dj(x) dΓ
                int = compute_face_integral_m2d_time_3dcoeffs(
                    mi,facedata.dΩ,coeffs3d,get_hf(facedata)
                )


                dmat[local_pos,i] += int/get_area(facedata.dΩ)
            end 
        end 
    end
    nothing
end

function _handle_volume_moments!(
    dmat::AbstractMatrix,
    base3d,
    mesh::Mesh{3,ET},
    vol_id::Int,
    ntl,
    hvol::Float64,
    vol_data::VolumeIntegralData,
    abs_volume::Float64
) where {K,ET<:ElType{K}}
 
    K ≤ 1 && return nothing

    mombase3d = get_base(BaseInfo{3,K-2,1}()).base
    moment_ids = mesh.int_coords_connect[4][vol_id]

 
    for (i,mi) in enumerate(base3d)
        for (j,mj) in enumerate(mombase3d)
        
            moment_id = moment_ids[j]
            local_pos = get_local_id(ntl,moment_id)
            val = compute_volume_integral_unshifted(mj*mi,vol_data,hvol)
            dmat[local_pos,i] += val/abs_volume
        end 
        
    end
    nothing

end

function create_volume_dmat(volume_id::Int,
    mesh::Mesh{3,ET}, 
    facedata_col::Dict{Int,<:FaceData{3,L}},
    vol_data::VolumeIntegralData,
    volume_node_mapping::NodeID2LocalID) where {K,L,ET<:ElType{K}}

    hvol = vol_data.hvol 
    bcvol = vol_data.vol_bc 
    abs_volume = vol_data.integrals[1]

 
    n_nodes = get_n_nodes(volume_node_mapping)
    base3d  = get_base(BaseInfo{3,K,1}())

    dmat    = FixedSizeMatrix{Float64}(undef,n_nodes,length(base3d))
    dmat   .= 0.0
 
    for (node_id,local_pos) in volume_node_mapping.map

        !_is_vertex_like_node(node_id,mesh) && continue

        dmat[local_pos,:] .= base_eval(base3d,mesh[node_id],bcvol,hvol)
    end


    _handle_face_moments!(
        dmat,mesh,bcvol,facedata_col,volume_node_mapping,volume_id,hvol
    )
    _handle_volume_moments!(
        dmat,base3d,mesh,volume_id,volume_node_mapping,hvol,vol_data,abs_volume
    )

    return dmat
end


function LinearAlgebra.Matrix{T}(::Any,::Static.StaticInt{N},n::Int) where {T,N}
    return Matrix{T}(undef,N,n)
end


function create_volume_vem_projectors(
    volume_id::Int,
    mesh::Mesh{3,ET}, 
    vol_data::VolumeIntegralData,
    facedata_col::Dict{Int,<:FaceData{3,L}},
    volume_node_mapping::NodeID2LocalID) where {K,L,ET<:ElType{K}}

    @assert K <= 2 "currently only implemented for K=1 and K=2"

    hvol        = vol_data.hvol
    abs_volume  = vol_data.integrals[1]
    bcvol       = vol_data.vol_bc

    base3d = get_base(BaseInfo{3,K,1}()).base
    dmat = create_volume_dmat(volume_id,mesh,facedata_col,vol_data,volume_node_mapping)

    bmat = create_volume_bmat(volume_id,mesh,bcvol,hvol,abs_volume,facedata_col,volume_node_mapping)
    gmat = static_matmul(bmat,dmat,Val((length(base3d),length(base3d))))


    proj_s = FixedSizeMatrix{Float64}(undef,size(gmat,1),size(bmat,2))
    proj   = FixedSizeMatrix{Float64}(undef,size(dmat,1),size(proj_s,2))
    
    Octavian.matmul!(proj_s,inv(gmat),bmat)
    Octavian.matmul!(proj,dmat,proj_s)
    return proj_s, proj
end

