using Ju3VEM
using Ju3VEM.VEMGeo: project_to_3d_flat
using Bumper
"""
    get_outward_normal(bcvol,face_data::FaceData{3,K,L}) where {K,L}

Returns the outward normal of the face w.r.t the volume center. 
!!!WARNING!!! Only works reliable for convex volumes.
"""
function get_outward_normal(bcvol,face_data::FaceData{3,K,L}) where {K,L}
    #TODO: move this function closer to geometric utils of the package
    plane = face_data.dΩ.plane
    normal = plane.n 

    face_bc3d = project_to_3d(get_bc(face_data),plane)
    face_to_vol_center = bcvol - face_bc3d 

    if dot(face_to_vol_center,normal) < 0
        return normal
    else
        return -normal
    end
end

function get_outward_normal(bcvol,face_data::FaceData{2,K,L}) where {K,L}
    return face_data.dΩ.plane.n
end


function precompute_base_projection_coeffs(::BaseInfo{3,O},
    dΩ::FaceIntegralData,bc_vol::SVector{3,Float64}) where {O}

    plane = dΩ.plane
    u     = plane.u
    v     = plane.v
    # resets offset from the precomputed integrals
    bc3d_face  = project_to_3d(get_bc(dΩ),plane)

    # offset of the volume center 
    bc_vol2d   = project_to_2d_abs(bc_vol,plane)
    flat_bcvol = project_to_3d_flat(bc_vol2d,plane)

    base = get_base(BaseInfo{3,O-1,1}()).base 

    coeff_list = MVector{length(base),SVector{length(base),Float64}}(undef)

    @no_escape begin 
        vec = @alloc(Float64,length(base))
        for (i,m) in enumerate(base)
            vec .= 0.0 
            compute_transformation_coeffs3d_to_2d!(
                vec,m,(bc3d_face-flat_bcvol),u,v,Val(O))

            coeff_list[i] = SVector{length(base),Float64}(vec)
        end
    end

    return SVector(coeff_list)
end

function get_dot_product_projection_coeffs(∇m::SVector{3,Monomial{Float64,3}}, 
    n::SVector{3,Float64},
    coeff_list::SVector{L,SVector{L,Float64}},hv) where {L} 


    coeffs = zero(eltype(coeff_list))

    for (i,∇mi) in enumerate(∇m)
        order = sum(∇mi.exp)
        idx = get_exp_to_idx_dict(∇mi.exp)
        coeffs += coeff_list[idx]*n[i]*∇mi.val/hv^order
    end

    return coeffs
end


function compute_face_integral_m2d_time_3dcoeffs(m2d::Monomial{T,2},
    fd::FaceIntegralData{D,K,L},coeffs::SVector,
    hf::T = fd.hf,hv::T = 1.0) where {D,T,K,L}


    order2 = sum(m2d.exp)

    full_base = get_base(BaseInfo{2,K,1}()).base


    ∫m = zero(T)
    for i in eachindex(coeffs)
        ci = coeffs[i]
        new_exp   = m2d.exp + full_base[i].exp 
        new_coeff = m2d.val*ci 

        idx = get_exp_to_idx_dict(new_exp)
        ∫m += new_coeff*fd.integrals[idx]
    end
    ∫m/(hf^order2)

    return ∫m
end



function create_volume_bmat(volume_id::Int,
    mesh::Mesh{3,ET},
    bcvol,hvol,
    facedata_col::Dict{Int,FaceData{3,K,L}},
    volume_node_mapping::NodeID2LocalID) where {K,L,ET<:ElType{K}}


    topo = mesh.topo

    base_2d = get_base(BaseInfo{2,K,1}()).base
    base_3d = get_base(BaseInfo{3,K,1}()).base


    n_nodes = get_n_nodes(volume_node_mapping)
    bmat = zeros(length(base_3d),n_nodes)

    if K == 1 
        bmat[1,:] .= 1/n_nodes
    end

    iterate_volume_areas(facedata_col,topo,volume_id) do _, face_data, _

        face_normal      = get_outward_normal(bcvol,face_data)
        ΠsL2             = face_data.ΠsL2
        face_node_ids    = face_data.face_node_ids

        hf = face_data.dΩ.hf

        coeff_list = precompute_base_projection_coeffs(BaseInfo{3,K,1}(),
                        face_data.dΩ,bcvol)

        integrated_vals = MVector{length(base_2d),Float64}(undef)
        for (i,m3d) in enumerate(base_3d)
            i == 1 && continue # gradient of constant is zero
            ∇m3d = ∇(m3d,hvol)
            integrated_vals .= 0.0
            for (j,m2d) in enumerate(base_2d)
                ∇m3d_dot_n = get_dot_product_projection_coeffs(
                    ∇m3d,face_normal,coeff_list,hvol)

                integrated_vals[j] += compute_face_integral_m2d_time_3dcoeffs(
                    m2d,face_data.dΩ,∇m3d_dot_n,hf)
            end


            for (nc,colv) in enumerate(eachcol(ΠsL2))
                val = dot(colv,integrated_vals)
                col_id  = get_local_id(volume_node_mapping,face_node_ids[nc])
                bmat[i,col_id] += val
            end
        end
    end
    return bmat
end



