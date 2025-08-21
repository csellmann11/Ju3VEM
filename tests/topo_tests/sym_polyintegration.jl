using Ju3VEM
using Ju3VEM.VEMGeo: get_base_len,compute_transformation_coeffs3d_to_2d!,get_plane_parameters
using Ju3VEM.VEMGeo: compute_transformation_coeffs2d_to_2d!,compute_transformation_coeffs2d_to_2d
using StaticArrays, LoopVectorization 
using LinearAlgebra, FastGaussQuadrature
using Symbolics, Chairmarks
using Test
using FixedSizeArrays, Bumper

const __GAUSS_LEGENDRE_PW = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()

function get_gauss_legendre(n::Int)
    get!(__GAUSS_LEGENDRE_PW,n) do
        FastGaussQuadrature.gausslegendre(n)
    end
end


struct FaceData{D,K,L} 
    hf::Float64
    bc::SVector{2,Float64} 
    p0::SVector{D,Float64}
    u ::SVector{D,Float64} 
    v ::SVector{D,Float64}  
    n ::SVector{D,Float64}

    integrals::SVector{L,Float64}
    # maybe add the projection matrix here?
end

function FaceData{K}(hf::Float64,bc::SVector{2,Float64},p0::SVector{D,Float64},
    u::SVector{D,Float64},v::SVector{D,Float64},
    n::SVector{D,Float64},integrals::SVector{L,Float64}) where {D,K,L}

    return FaceData{D,K,L}(hf,bc,p0,u,v,n,integrals)
end


function face2d_sym_integral(
    node_coords::AbstractVector{<:StaticVector{D,T}}, 
    e1::Int,e2::Int,
    p0::StaticVector{D,T} = SA[0.0,0.0],
    u::StaticVector{D,T} = SA[1.0,0.0],
    v::StaticVector{D,T} = SA[0.0,1.0]) where {D,T}

    int_val = zero(T)

    ngpoints = cld(e1 + e2 + 2, 2)
    gpoints,gweights =  
            get_gauss_legendre(ngpoints)


    for i in eachindex(node_coords)
        ip1 = get_next_idx(node_coords,i)
        
        nci2d   = SA[dot(node_coords[i]-p0,u)  ,dot(node_coords[i]-p0,v)] 
        nci2d_n = SA[dot(node_coords[ip1]-p0,u),dot(node_coords[ip1]-p0,v)]

        x1,y1 = nci2d   
        x2,y2 = nci2d_n 


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
        
    end
    int_val /= (e1+1)
    return -int_val
end



"""
    Function precomputes the integrals of the monomials on a face up to the monomial order K
"""
function precompute_face_monomials(face_nodes::AbstractVector{<:StaticVector{D,T}},::Val{_K}) where {D,T,_K}

    K    = max(_K,1)
    base = get_base(BaseInfo{2,K,1}()).base  

    hf           = max_node_distance(face_nodes)
    u,v,n,_,p0   = get_plane_parameters(face_nodes) 
    p02d         = SA[dot(p0,u),dot(p0,v)]

    L         = length(base)
    integrals = MVector{L,T}(undef)

    for i in 1:L
        m   = base[i]
        exp = m.exp
        integrals[i] = face2d_sym_integral(face_nodes,exp[1],exp[2],p0,u,v)
    end

    face_area = integrals[1] 
    bc        = SA[integrals[2],integrals[3]]/face_area + p02d
    return FaceData{K}(hf,bc,p0,u,v,n,SVector(integrals))
end

function precompute_face_monomials(face_id::Int,topo::Topology{D},::Val{K}) where {D,K}
    face_node_ids = get_area_node_ids(topo,face_id)
    face_nodes    = @view topo.nodes[face_node_ids] 
    return precompute_face_monomials(face_nodes,Val(K))
end

function compute_face_integral(m::Monomial{T,3},fd::FaceData{D,K,L},
    bc3d::SVector{D,T} = fd.u*fd.bc[1] + fd.v*fd.bc[2],
    h::T = fd.hf) where {D,T,K,L}

    order = sum(m.exp) 
    len   = get_base_len(2,order,1)

    @assert len ≤ length(fd.integrals) "not enough integrals where precomputed"

    ∫m = @no_escape begin 
        coeffs = @alloc(Float64,len)
        compute_transformation_coeffs3d_to_2d!(
            coeffs,m, fd.p0-bc3d,fd.u,fd.v)


        ∫m = zero(T)
        @inbounds for (i,ci) in enumerate(coeffs)
            ∫m += ci*fd.integrals[i]
        end
        ∫m/(h^order)
    end
    return ∫m
end


function compute_face_integral(m::Monomial{T,2},fd::FaceData{D,K,L},
    bc2d::SVector{2,T} = fd.bc,
    h::T = fd.hf) where {D,T,K,L}

    order = sum(m.exp) 
    len   = get_base_len(2,order,1)
    p02d  = SA[dot(fd.p0,fd.u),dot(fd.p0,fd.v)]

    @assert len ≤ length(fd.integrals) "not enough integrals where precomputed"

    ∫m = @no_escape begin 
        coeffs = @alloc(Float64,len)
        compute_transformation_coeffs2d_to_2d!(coeffs,m,p02d-bc2d)

        ∫m = zero(T)
        @inbounds for (i,ci) in enumerate(coeffs)
            ∫m += ci*fd.integrals[i]
        end
        ∫m/(h^order)
    end
    return ∫m
end







quad_nodes = [SA[0.0,0.0,0.0], SA[1.0,0.0,0.0], SA[1.0,1.0,1.0], SA[0.0,1.0,1.0]] .* 2 .+ Scalar(SA[1.0,1.0,1.0])
hf = max_node_distance(quad_nodes)
u, v, n, nus,p0 = get_plane_parameters(quad_nodes)

p02d = SA[dot(p0,u),dot(p0,v)]

a   = face2d_sym_integral(quad_nodes,0,0,p0,u,v)
bcu = face2d_sym_integral(quad_nodes,1,0,p0,u,v)/a
bcv = face2d_sym_integral(quad_nodes,0,1,p0,u,v)/a

face_data = precompute_face_monomials(quad_nodes,Val(5))

face_data.integrals[1]/a

@test face_data.bc ≈ SA[bcu,bcv] + p02d


x1 = quad_nodes[1]; x2 = quad_nodes[2]; x3 = quad_nodes[3]; x4 = quad_nodes[4]
x2d1 = SA[dot(x1-p0,u),dot(x1-p0,v)]
x2d2 = SA[dot(x2-p0,u),dot(x2-p0,v)]
x2d3 = SA[dot(x3-p0,u),dot(x3-p0,v)]
x2d4 = SA[dot(x4-p0,u),dot(x4-p0,v)]

quad_nodes_2d = [x2d1, x2d2, x2d3, x2d4] 

a2d   = face2d_sym_integral(quad_nodes_2d,0,0)
bcu2d = face2d_sym_integral(quad_nodes_2d,1,0)/a
bcv2d = face2d_sym_integral(quad_nodes_2d,0,1)/a
bc    = SA[bcu,bcv] 
bc3d = SA[2.0,2.0,2.0]


@test bcu ≈ bcu2d 
# @test bcv ≈ bcv2d + p02d[2] #+ p02d[2]

bc3d = (face_data.bc[1]-p02d[1])*u + (face_data.bc[2]-p02d[2])*v + p0
@test bc3d ≈ SA[2.0,2.0,2.0]

@test x2 ≈ x2d2[1]*u + x2d2[2]*v + x1
@test x3 ≈ x2d3[1]*u + x2d3[2]*v + x1
@test x4 ≈ x2d4[1]*u + x2d4[2]*v + x1

true_area = (x2d2[1] - x2d1[1])*(x2d3[2] - x2d1[2])


m3d = Monomial(1.0,SA[2,1,1])
coeffs2d = compute_transformation_coeffs3d_to_2d(m3d,p0,u,v)
poly2d   = Polynomial(coeffs2d,get_base(BaseInfo{2,4,1}()).base)


bc3d = face_data.u*(face_data.bc[1] - p02d[1]) + face_data.v*(face_data.bc[2] - p02d[2]) + face_data.p0




I = compute_face_integral(m3d,face_data)

ref = 4*sqrt(2)/(9*hf^4)

@test I ≈ ref






quad_nodes2 = [SA[0.0,0.0,0.0], SA[1.0,0.0,0.0], SA[1.0,1.0,0.0], SA[0.0,1.0,0.0]] .* 2 .+ Scalar(SA[1.0,1.0,1.0])
hf = max_node_distance(quad_nodes2)
u, v, n, nus,p0 = Ju3VEM.VEMGeo.get_plane_parameters(quad_nodes2)

# p0 = 0*p0
x1d2 = SA[dot(quad_nodes2[1],u),dot(quad_nodes2[1],v)]
x2d2 = SA[dot(quad_nodes2[2],u),dot(quad_nodes2[2],v)]
x2d3 = SA[dot(quad_nodes2[3],u),dot(quad_nodes2[3],v)]
x2d4 = SA[dot(quad_nodes2[4],u),dot(quad_nodes2[4],v)]


_x1d2 = SA[dot(quad_nodes2[1]-p0,u),dot(quad_nodes2[1]-p0,v)]
_x2d2 = SA[dot(quad_nodes2[2]-p0,u),dot(quad_nodes2[2]-p0,v)]
_x2d3 = SA[dot(quad_nodes2[3]-p0,u),dot(quad_nodes2[3]-p0,v)]
_x2d4 = SA[dot(quad_nodes2[4]-p0,u),dot(quad_nodes2[4]-p0,v)]

p02d = SA[dot(p0,u),dot(p0,v)]

@test x1d2 == _x1d2 + p02d
@test x2d2 == _x2d2 + p02d
@test x2d3 == _x2d3 + p02d
@test x2d4 == _x2d4 + p02d




face_data2 = precompute_face_monomials(quad_nodes2,Val(5))
bc3d = face_data2.u*(face_data2.bc[1] - p02d[1]) + face_data2.v*(face_data2.bc[2] - p02d[2]) + face_data2.p0

m3d = Monomial(1.0,SA[0,1,0])
compute_transformation_coeffs3d_to_2d(m3d,p0,u,v)
I = compute_face_integral(m3d,face_data2,zero(bc3d),2.0)



face_data2.integrals[1:3]

m2d = Monomial(1.0,SA[0,1])

compute_transformation_coeffs2d_to_2d(m2d,p02d)
I2 = compute_face_integral(m2d,face_data2,zero(p02d),2.0)

@test I ≈ I2









quad_nodes3 = [SA[0.0,0.0,0.0], SA[1.0,0.0,0.0], SA[1.0,1.0,0.0], SA[0.0,1.0,0.0]] .* 2 .+ 2*Scalar(SA[1.0,1.0,0.0])
face_data3 = precompute_face_monomials(quad_nodes3,Val(5))

m3d = Monomial(1.0,SA[0,1,0])
m2d = Monomial(1.0,SA[0,1])

I  = compute_face_integral(m3d,face_data3,zero(bc3d),2.0)
I2 = compute_face_integral(m2d,face_data3,zero(p02d),2.0)

b = @b compute_face_integral($m2d,$face_data3,zero($p02d),2.0)
display(b)

# @test I  ≈ 6.0 
# @test I2 ≈ 6.0 








