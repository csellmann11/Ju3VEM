# using Ju3VEM
# using Ju3VEM.VEMGeo: get_base_len,compute_transformation_coeffs3d_to_2d!,get_plane_parameters
# using Ju3VEM.VEMGeo: compute_transformation_coeffs2d_to_2d!,compute_transformation_coeffs2d_to_2d
# using StaticArrays, LoopVectorization 
# using LinearAlgebra, FastGaussQuadrature
# using Symbolics, Chairmarks
# using Test
# using FixedSizeArrays, Bumper

const __GAUSS_LEGENDRE_PW = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()
const __GAUSS_LOBATTO_PW  = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()

function get_gauss_legendre(n::Int)
    get!(__GAUSS_LEGENDRE_PW,n) do
        FastGaussQuadrature.gausslegendre(n)
    end
end

function get_gauss_lobatto(n::Int)
    get!(__GAUSS_LOBATTO_PW,n) do
        FastGaussQuadrature.gausslobatto(n)
    end
end

struct D2FaceParametrization{D}
    p0      ::SVector{D,Float64}
    p0_2d   ::SVector{2,Float64}
    u       ::SVector{D,Float64} 
    v       ::SVector{D,Float64}  
    n       ::SVector{3,Float64}
end

function D2FaceParametrization(
    face_nodes::AbstractVector{<:StaticVector{3,T}}) where {T}

    u,v,n,_,p0   = get_plane_parameters(face_nodes) 
    p0_2d = SA[dot(p0,u),dot(p0,v)]
    return D2FaceParametrization(SVector(p0),p0_2d,u,v,n)
end

function D2FaceParametrization(
    x::AbstractVector{<:StaticVector{2,T}}) where {T}

    u = SA[1.0,0.0]
    v = SA[0.0,1.0]
    n = SA[0.0,0.0,1.0]
    p0 = SVector(x[1])

    return D2FaceParametrization(p0,p0,u,v,n)
end

function project_to_2d_abs(
    x::StaticVector{3,T},plane::D2FaceParametrization{D}) where {T,D}

    x2d = SA[dot(x,plane.u),dot(x,plane.v)] 
    return x2d
end

@inline project_to_2d_abs(x::StaticVector{2},
    ::D2FaceParametrization{2}) = x

function project_to_2d_rel(
    x::StaticVector{3,T},plane::D2FaceParametrization{D}) where {T,D}

    x2d = SA[dot(x,plane.u),dot(x,plane.v)] - plane.p0_2d
    return x2d
end

@inline project_to_2d_rel(x::StaticVector{2},
    plane::D2FaceParametrization{2}) = x - plane.p0_2d


function project_to_3d(
    x::StaticVector{2,T},plane::D2FaceParametrization{3}) where {T}

    x3d = plane.u*(x[1] - plane.p0_2d[1]) + 
          plane.v*(x[2] - plane.p0_2d[2]) + plane.p0
    return x3d
end




"""
    project_to_3d_flat(x::StaticVector{2,T},plane::D2FaceParametrization{3}) where {T}

    Function does in almost all cases the same as project_to_3d. 
    It only differs if the normal of the plane is parallel to 
    on of the coordinate axes.
"""
function project_to_3d_flat(
    x::StaticVector{2,T},plane::D2FaceParametrization{3}) where {T}

    x3d = plane.u*(x[1]) + 
          plane.v*(x[2]) 
    return x3d
end




struct FaceIntegralData{D,L}  
    hf::Float64 
    bc::SVector{2,Float64} 
    plane::D2FaceParametrization{D}
 
    integrals    ::SVector{L,Float64}
end 
 
# function FaceIntegralData(hf::Float64,bc::SVector{2,Float64},
#     plane::D2FaceParametrization{D},integrals::SVector{L,Float64}) where {D,L}


#     return FaceIntegralData{D,L}(hf,bc,plane,integrals)
# end

@inline get_area(fd::FaceIntegralData) = fd.integrals[1]
@inline get_bc(fd::FaceIntegralData)   = fd.bc


struct FaceData{D,L,FV<:FlattenVecs{3,Int}} 
    face_node_ids::FV 
    dΩ           ::FaceIntegralData{D,L}
    ΠsL2         ::FixedSizeMatrixDefault{Float64}
end 
@inline get_area(fd::FaceData) = get_area(fd.dΩ)
@inline get_bc(fd::FaceData)   = get_bc(fd.dΩ)
@inline get_hf(fd::FaceData)   = fd.dΩ.hf

"""
    get_outward_normal(bcvol,face_data::FaceData{3,L}) where {L}

Returns the outward normal of the face w.r.t the volume center. 
!!!WARNING!!! Only works reliable for convex volumes.
"""
function get_outward_normal(bcvol,face_data::FaceData{3,L}) where {L}
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


function face2d_sym_integral(
    node_coords::AbstractVector{<:StaticVector{D,T}}, 
    e1::Int,e2::Int,bc::StaticVector{2,T},
    plane::D2FaceParametrization{D}) where {D,T}

    int_val = zero(T) 

    ngpoints = cld(e1 + e2 + 2, 2)
    gpoints,gweights =  
            get_gauss_legendre(ngpoints)

    nci2d   = project_to_2d_abs(node_coords[1],plane) 
    for i in eachindex(node_coords)
        ip1 = get_next_idx(node_coords,i)
        
        nci2d_n = project_to_2d_abs(node_coords[ip1],plane)

        x1,y1 = nci2d   - bc
        x2,y2 = nci2d_n - bc


        e1_int = zero(T)
        @turbo for j in eachindex(gpoints,gweights) # calc int x^(e1+1)y^e2 dxi
            ξ  = gpoints[j]
            w  = gweights[j]

            f1 = 0.5(1-ξ)
            f2 = 0.5(1+ξ) 


            x_intp = (f1*x1+f2*x2)
            y_intp = (f1*y1+f2*y2)

            val = x_intp^(e1+1)*y_intp^e2*w
            e1_int += val
        end
        # (y1-y2) = nx * L
        int_val += (y1-y2)*e1_int/2
        nci2d   = nci2d_n
    end
    int_val /= (e1+1)
    return -int_val
end



"""
    Function precomputes the integrals of the monomials on a face up to the monomial order K
"""
function precompute_face_monomials(
    face_nodes::AbstractVector{<:StaticVector{D,T}},::Val{K_MAX}) where {D,T,K_MAX}

    K    = max(K_MAX,1)
    base = get_base(BaseInfo{2,K,1}()).base  

    hf           = max_node_distance(face_nodes)
    plane = D2FaceParametrization(face_nodes)


    L         = length(base)
    integrals = zero(MVector{L,T})
    face_area = face2d_sym_integral(face_nodes,0,0,SA[0.0,0.0],plane)
    bcu       = face2d_sym_integral(face_nodes,1,0,SA[0.0,0.0],plane)/face_area
    bcv       = face2d_sym_integral(face_nodes,0,1,SA[0.0,0.0],plane)/face_area
    bc        = SA[bcu,bcv]
    integrals[1] = face_area

    for i in 4:L
        m   = base[i]
        exp = m.exp
        integrals[i] = face2d_sym_integral(face_nodes,exp[1],exp[2],bc,plane)
    end

    return FaceIntegralData(hf,bc,plane,SVector(integrals))
end 

function precompute_face_monomials(
    face_id::Int,topo::Topology{D},::Val{K}) where {D,K}

    face_node_ids = get_area_node_ids(topo,face_id)
    face_nodes    = @view topo.nodes[face_node_ids] 
    return precompute_face_monomials(face_nodes,Val(K))
end


function compute_face_integral_coeffs!(
    coeffs::AbstractVector{<:Float64},
    m::Monomial{T,3},
    fd::FaceIntegralData{D,L},
    offset::BT = nothing
    ) where {D,T,L,BT<:Union{Nothing,SVector{3,T}}}

    plane = fd.plane 
 
    bc3d_face  = project_to_3d(fd.bc,plane)

    offset_2d = if offset === nothing 
        project_to_2d_abs(bc3d_face,plane)
    else 
        project_to_2d_abs(offset,plane)
    end

    flat_offset = project_to_3d_flat(offset_2d,plane)


    compute_transformation_coeffs3d_to_2d!(
        coeffs,m, (bc3d_face-flat_offset),plane.u,plane.v)

    return coeffs
end


function compute_face_integral(m::Monomial{T,3},
    fd::FaceIntegralData{D,L},
    offset::BT = nothing,
    h::T = fd.hf) where {D,T,L,BT<:Union{Nothing,SVector{3,T}}}

    order = sum(m.exp) 
    len   = get_base_len(2,order,1)

    @assert len ≤ length(fd.integrals) "not enough integrals where precomputed"

    ∫m = @no_escape begin 
        coeffs = @alloc(Float64,len)
        compute_face_integral_coeffs!(coeffs,m,fd,offset)

        ∫m = zero(T)
        @inbounds for (i,ci) in enumerate(coeffs)
            ∫m += ci*fd.integrals[i]
        end
        ∫m/(h^order)
    end
    return ∫m
end


function compute_face_integral(m::Monomial{T,2},fd::FaceIntegralData{D,L},
    offset2d::SVector{2,T} = fd.bc,
    h::T = fd.hf) where {D,T,L}

    order = sum(m.exp) 
    len   = get_base_len(2,order,1)

    plane = fd.plane

    @assert len ≤ length(fd.integrals) "not enough integrals where precomputed"

    ∫m = @no_escape begin 
        coeffs = @alloc(Float64,len)
        compute_transformation_coeffs2d_to_2d!(coeffs,m,fd.bc-offset2d)

        ∫m = zero(T)
        @inbounds for (i,ci) in enumerate(coeffs)
            ∫m += ci*fd.integrals[i]
        end
        ∫m/(h^order)
    end
    return ∫m
end

function compute_face_integral_unshifted(m::Monomial{T,2},
    fd::FaceIntegralData{D,L},h::T = fd.hf) where {D,T,L}
    idx   = get_exp_to_idx_dict(m.exp)
    order = sum(m.exp)
    return fd.integrals[idx]/(h^order)
end



# ================================================00
# Volume integration 

struct VolumeIntegralData{L}
    hvol::Float64
    vol_bc::SVector{3,Float64}
    integrals::SVector{L,Float64}
end 


"""
    precompute_volume_monomials(vol_id::Int,topo::Topology,facedata_col::Dict{Int,<:FaceData{3}},::Val{K}) where {K}

Precomputes the integrals of the monomials on a volume up to the monomial order K.


# WARNING: currently requires the volume to be convex (_some_point_inside) must be actually inside the volume

"""
function precompute_volume_monomials(vol_id::Int,
    topo::Topology,
    facedata_col::Dict{Int,<:FaceData{3}},::Val{K_MAX}, 
    shift_bc = true) where {K_MAX}



    base = get_base(BaseInfo{3,K_MAX,1}()).base

    L = length(base)
    integrals = zero(MVector{L,Float64})

    hvol = max_node_distance(topo.nodes,get_volume_node_ids(topo,vol_id))

    # Maby replace this with the mean of the face barycenters
    _some_point_inside = mean(get_nodes(topo)[ni] for ni in get_volume_node_ids(topo,vol_id))

    geo_data = zero(MVector{4,Float64}) 
    # the areea and th ebary cneter are computed without the offset 
    # therefore they are in an extra loop

    #TODO: swith loop order for performance
    for i in 1:min(4,L)
        mi = base[i] 
        p = mi.exp[1]
 
        m_face = Monomial(1/(p+1),SA[p+1,mi.exp[2],mi.exp[3]])

        iterate_volume_areas(facedata_col,topo,vol_id) do _, facedata, _
            dΩ = facedata.dΩ 
            normal = get_outward_normal(_some_point_inside,facedata)
            geo_data[i] += compute_face_integral(m_face*normal[1],dΩ,SA[0.0,0.0,0.0],1.0)
        end
    end
   
    volume_measure = geo_data[1]
    vol_bc         = SA[geo_data[2],geo_data[3],geo_data[4]]/volume_measure

    integrals[1] = volume_measure
    if !shift_bc 
        integrals[2:4] .= SA[geo_data[2],geo_data[3],geo_data[4]]
    end

    integral_center = shift_bc ? vol_bc : zero(vol_bc)

    for i in 5:length(base)
        mi = base[i]
        p = mi.exp[1]
        m_face = Monomial(1/(p+1),SA[p+1,mi.exp[2],mi.exp[3]])


        
        iterate_volume_areas(facedata_col,topo,vol_id) do _, facedata, _
            dΩ = facedata.dΩ 
            normal = get_outward_normal(_some_point_inside,facedata)
            integrals[i] += compute_face_integral(m_face*normal[1],dΩ,integral_center,1.0)#hvol)
        end
    end

    return VolumeIntegralData(hvol,vol_bc,SVector(integrals))
end


function compute_volume_integral_unshifted(m::Monomial{T,3},
    vol_data::VolumeIntegralData{L},h::T = vol_data.hvol) where {T,L}

    idx   = get_exp_to_idx_dict(m.exp)
    order = sum(m.exp)
    return m.val*vol_data.integrals[idx]/(h^order)
end

















