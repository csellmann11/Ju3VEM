import Ferrite as FR

function build_kel_quad_point!(
    kelement::CachedMatrix{Float64},
    rhs_element::CachedVector{Float64},
    cv::CellValues{D,U,ET},
    fe_cv::FR.CellValues,
    Nv::AbstractVector,
    ∇Nv::AbstractVector,
    f::F,
    qpoint::Integer,
    x::StaticVector{3,Float64},
    mat_law::Helmholtz,
    pars...) where {D,U,ET,F<:Function}

 
    hvol = cv.volume_data.hvol
    bc = cv.volume_data.vol_bc
    x_scaled = (x .- bc)/hvol

    dΩ = FR.getdetJdV(fe_cv,qpoint) 

    #INFO: for linea materials the input grad is just a dummy
    ∇Nvx = [∇Nv[i](x_scaled) for i in eachindex(Nv)]
    ℂ = eval_hessian(mat_law,∇Nvx[1],pars...)
    



    # @inbounds for (i,∇Ni) in enumerate(∇Nv)
    #     Nix  = Nv[i](x_scaled)
    #     ∇Nixℂ = ℂ ⊡₂ ∇Ni(x_scaled) 
    #     rhs_element[i] += f(x)⋅Nix*dΩ
    #     for j in i:length(Nv)
    #         ∇Njx = ∇Nv[j](x_scaled)
    #         kelement[i,j] += ∇Nixℂ ⋅ ∇Njx *dΩ
    #     end
    # end
    @inbounds for (i,∇Nix) in enumerate(∇Nvx)
        rhs_element[i] +=  f(x)⋅Nv[i](x_scaled)*dΩ
        ∇Nixℂ = ℂ ⊡₂ ∇Nix
        for j in i:length(Nv)
            kelement[i,j] += ∇Nixℂ ⋅ ∇Nvx[j] *dΩ
        end
    end


    kelement
end

function symmetrize_kel!(
    kelement::CachedMatrix{Float64})

    @assert size(kelement,1) == size(kelement,2) "kelement must be square"
    @inbounds for i in axes(kelement,1)
        for j in (i+1):size(kelement,2)  
            kelement[j,i] = kelement[i,j]  
        end
    end
    kelement
end



function tetrahedron_volume(nodes::SVector{4, <:FR.Vec})
    # Volume = |det([p1-p4, p2-p4, p3-p4])| / 6
    v1 = nodes[1] - nodes[4]
    v2 = nodes[2] - nodes[4]
    v3 = nodes[3] - nodes[4]
    
    # Compute the scalar triple product (determinant)
    return abs(v1 ⋅ (v2 × v3)) / 6
end
function build_local_kel_and_f!(
    kelement::CachedMatrix{Float64},
    rhs_element::CachedVector{Float64},
    cv::CellValues{D,U,ET},
    element_id::Integer,
    f::F,
    fe_cv::FR.CellValues,
    mat_law::Helmholtz,
    γ::Float64 = 1/4,
    pars...) where {D,U,F<:Function,K,ET<:ElType{K}}


    hvol = cv.volume_data.hvol

    base3d = get_base(BaseInfo{3,K,U}())

    proj_s, proj = create_volume_vem_projectors(
        element_id,cv.mesh,cv.volume_data,cv.facedata_col,cv.vnm)

    ebf = ElementBaseFunctions(base3d,stretch(proj_s,Val(U)))

    poly_type = eltype(ebf)
    grad_type = Core.Compiler.return_type(∇p, Tuple{poly_type,Float64})

  
    ∇Nv = FixedSizeVector{grad_type}(undef,length(ebf))
    Nv  = FixedSizeVector{poly_type}(undef,length(ebf))

    for (i,p) in enumerate(ebf)
        ∇Nv[i] = ∇p(p,hvol)
        Nv[i]  = p |> poly_type
    end

    tets_local, l2g = tetrahedralize_volume(cv.mesh.topo, element_id)

    setsize!(kelement,(length(ebf),length(ebf)))
    setsize!(rhs_element,(length(ebf),))

    ref_volume = cv.volume_data.integrals[1]
    test_volume = 0.0

    for tet_local_ids in tets_local
        tets_global_ids = SVector{4,Int32}(l2g[id] for id in tet_local_ids)
        tet_nodes = FR.Vec{D}.(get_nodes(cv.mesh.topo)[tets_global_ids])
        FR.reinit!(fe_cv,tet_nodes)

        for qpoint in 1:FR.getnquadpoints(fe_cv)
            
            x = FR.spatial_coordinate(fe_cv,qpoint,tet_nodes) |> SVector{3,Float64}
            build_kel_quad_point!(kelement,rhs_element,cv,
                        fe_cv,Nv,∇Nv,f,qpoint,x,mat_law,pars...)

            test_volume += FR.getdetJdV(fe_cv,qpoint)
        end 
    end
    symmetrize_kel!(kelement)

    stab       = (I-proj)'*(I-proj)*hvol*γ
    kelement .+= stretch(stab,Val(U))

    return proj_s, proj
end





function build_local_kel_and_f!(
    kelement::CachedMatrix{Float64},
    rhs_element::CachedVector{Float64},
    cv::CellValues{D,U,ET},
    element_id::Integer,
    f::F,
    mat_law::Helmholtz,
    γ::Float64 = 1/4,
    pars::P = ()) where {D,U,F<:Function,ET<:ElType{1},P}


    hvol = cv.volume_data.hvol
    bc   = cv.volume_data.vol_bc
    dΩ   = cv.volume_data.integrals[1]

    base3d = get_base(BaseInfo{3,1,U}())
    L = length(base3d)


    k_poly_space = MMatrix{L,L,Float64} |> zero

    
    proj_s, proj = create_volume_vem_projectors(
        element_id,cv.mesh,cv.volume_data,cv.facedata_col,cv.vnm)

    n_dofs = size(proj,1)*U

    setsize!(kelement,(n_dofs,n_dofs))
    setsize!(rhs_element,(n_dofs,))


    ℂ0 = eval_hessian(mat_law,zero(SMatrix{U,D,Float64,U*D}),pars)
    ∇x_base = SVector{length(base3d)}(∇(m,hvol)(zero(bc)) for m in base3d)



    for (i,∇mx) in enumerate(∇x_base)
        ∇mxℂ0 = ℂ0 ⊡₂ ∇mx
        for (j,∇nx) in enumerate(∇x_base)
            k_poly_space[i,j] += ∇mxℂ0 ⋅ ∇nx * dΩ
        end
    end

    # computes proj_s' * k_poly_space * proj_s
    _temp = k_poly_space * stretch(proj_s,Val(U))
    mul!(kelement,stretch(proj_s',Val(U)),_temp)


    fmean = f(bc)
    for i in size(proj,1)
        phi_val = mean(proj[i,j] for j in axes(proj,2))
        for u in 1:U
            idx = (i-1)*U + u
            rhs_element[idx] += fmean[u] * phi_val * dΩ
        end
    end


    _Δmat = (I-proj)
    stab       = _Δmat'*_Δmat*hvol*γ
    kelement .+= stretch(stab,Val(U))

    return proj_s, proj
end
