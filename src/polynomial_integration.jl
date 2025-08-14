using FastGaussQuadrature
using StaticArrays
using LinearAlgebra
using LoopVectorization

const __GAUSS_LEGENDRE_PW = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()

function get_gauss_legendre(n::Int)
    get!(__GAUSS_LEGENDRE_PW,n) do
        FastGaussQuadrature.gausslegendre(n)
    end
end



"""
    steger_integral(node_coords::AbstractVector{V1}, e1::Int,e2::Int,bc::V2) where {V1<:AbstractVector{<:Real},V2<:AbstractVector{<:Real}}

Computes the integral of the polynomial x^e1*y^e2 over the plane defined by `node_coords`.
Uses the divergence theorem to compute the integral over the face. 
∫_A ∇⋅(x^(e1+1)*y^e2,0)*1/(e1+1)*dA = ∫_A x^e1*y^e2*dA = ∮_E 1/(e1+1)*x^(e1+1)*y^e2*nx*ds
"""
function face_integral(
    node_coords::AbstractVector{V1}, 
    e1::Int,e2::Int,bc::V2 = SVector(0.0,0.0)) where {T<:Real,V1<:AbstractVector{T},V2<:AbstractVector{T}}

    D = length(V1)
    @assert D == length(V2)

 
    int_val =zero(T)

    ngpoints = ceil(Int,(e1+e2+2)/2)
    gpoints,gweights =  
            get_gauss_legendre(ngpoints)

    for i in eachindex(node_coords)

        ip1 = get_next_idx(node_coords,i)
        
        x1,y1 = node_coords[i] - bc
        x2,y2 = node_coords[ip1]-bc


        e1_int = zero(T)
        for j in eachindex(gpoints,gweights) # calc int x^(e1+1)y^e2 dxi
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

function integral(m::Monomial{T,2},
    node_coords::AbstractVector,h::Real=1.,bc::AbstractVector = SVector(0.0,0.0)) where {T}

    e1,e2 = m.exp
    return m.val*face_integral(node_coords,e1,e2,bc)/h^(e1+e2)
end


# # ########## 3D Face Integration (planar polygon in 3D) ##########

# # Access Triangulation3D.project_to_plane without imposing strict include order
# try
#     using Main.Triangulation3D: project_to_plane
# catch
#     try
#         import ..Triangulation3D: project_to_plane
#     catch
#         # Fallback: local minimal implementation (kept in sync with src/triangulation.jl)
#         # Returns (origin, u, v, projected_points)
#         function project_to_plane(points::AbstractVector{<:SVector{3,T}}) where T
#             @assert length(points) ≥ 3
#             o = points[1]
#             # find non-identical point
#             i2 = 2
#             while i2 ≤ length(points) && isapprox(norm(points[i2] - o), zero(T); atol=eps(T)*10)
#                 i2 += 1
#             end
#             @assert i2 ≤ length(points)
#             e1 = points[i2] - o
#             # find non-collinear point
#             n = zero(SVector{3,T})
#             i3 = i2 + 1
#             while i3 ≤ length(points)
#                 e2 = points[i3] - o
#                 n = cross(e1, e2)
#                 if norm(n) > sqrt(eps(T))
#                     break
#                 end
#                 i3 += 1
#             end
#             @assert i3 ≤ length(points)
#             n = n / norm(n)
#             u = e1 - dot(e1, n) * n
#             u = u / norm(u)
#             v = cross(n, u)
#             pts2 = Vector{SVector{2,T}}(undef, length(points))
#             for (i,p) in enumerate(points)
#                 d = p - o
#                 pts2[i] = SVector(dot(d, u), dot(d, v))
#             end
#             return o, u, v, pts2
#         end
#     end
# end

# @inline function _shoelace_signed_area2(points::AbstractVector{<:SVector{2,T}}) where T
#     s = zero(T)
#     @inbounds for i in eachindex(points)
#         j = i == length(points) ? 1 : i + 1
#         s += points[i][1]*points[j][2] - points[j][1]*points[i][2]
#     end
#     return s/2
# end

# # Expand (a0 + a1 s + a2 t)^e into a 2D polynomial in s,t
# function _expand_linear_power(a0::T, a1::T, a2::T, e::Int) where T
#     degs = e
#     coeffs = zeros(T, degs+1, degs+1) # coeffs[k+1,l+1] for s^k t^l
#     @inbounds for k in 0:degs
#         for l in 0:(degs-k)
#             # multinomial coefficient: e!/(k! l! (e-k-l)!) = binomial(e,k)*binomial(e-k,l)
#             c = binomial(degs, k) * binomial(degs-k, l)
#             coeffs[k+1, l+1] = c * (a1^k) * (a2^l) * (a0^(degs-k-l))
#         end
#     end
#     return coeffs
# end

# # Multiply two 2D polynomials represented as coefficient matrices
# function _poly2d_mul(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where T
#     sA, tA = size(A)
#     sB, tB = size(B)
#     C = zeros(T, sA + sB - 1, tA + tB - 1)
#     for i in 1:sA
#         for j in 1:tA
#             aij = A[i,j]
#             aij == 0 && continue
#             @simd for k in 1:sB
#                 for l in 1:tB
#                     C[i+k-1, j+l-1] += aij * B[k,l]
#                 end
#             end
#         end
#     end
#     return C
# end

# """
#     face_integral(node_coords::AbstractVector{<:SVector{3,<:Real}}, ex::Int, ey::Int, ez::Int,
#                   bc::SVector{3,<:Real} = SVector(0.0,0.0,0.0)) -> Real

# Compute ∫_A (x-bcx)^ex (y-bcy)^ey (z-bcz)^ez dA over a planar polygon `node_coords` in 3D.
# Uses projection to an orthonormal in-plane basis and reduces to 2D monomial integrals.
# """
# function face_integral(node_coords::AbstractVector{V1}, ex::Int, ey::Int, ez::Int,
#     bc::V2 = SVector(0.0,0.0,0.0)) where {T<:Real,V1<:SVector{3,T}, V2<:SVector{3,T}}

#     @assert ex ≥ 0 && ey ≥ 0 && ez ≥ 0


#     # Project to plane: p = o + s*u + t*v, u,v orthonormal
#     pts3 = SVector{3,T}.(node_coords)
#     o, u, v, pts2 = project_to_plane(pts3)

#     # Ensure CCW orientation of pts2 for robust positive 2D area integration
#     if _shoelace_signed_area2(pts2) < 0
#         pts2 = reverse(pts2)
#     end

#     # Build polynomial P(s,t) = Π_q (α0_q + α1_q s + α2_q t)^{e_q}
#     α0x, α1x, α2x = (o[1]-bc[1]), u[1], v[1]
#     α0y, α1y, α2y = (o[2]-bc[2]), u[2], v[2]
#     α0z, α1z, α2z = (o[3]-bc[3]), u[3], v[3]


#     Px = _expand_linear_power(T(α0x), T(α1x), T(α2x), ex)
#     Py = _expand_linear_power(T(α0y), T(α1y), T(α2y), ey)
#     Pz = _expand_linear_power(T(α0z), T(α1z), T(α2z), ez)

#     Pxy = _poly2d_mul(Px, Py)
#     P = _poly2d_mul(Pxy, Pz)

#     # Integrate term-by-term: ∫ s^i t^j dA over pts2 using 2D face_integral
#     sdeg, tdeg = size(P)
#     total = zero(T)
#     @inbounds for i in 0:(sdeg-1)
#         for j in 0:(tdeg-1)
#             c = P[i+1, j+1]
#             c == 0 && continue
#             total += c * face_integral(pts2, i, j, SVector{2,T}(0,0))
#         end
#     end
#     return total
# end

# function integral(m::Monomial{T,3}, node_coords::AbstractVector{<:SVector{3}},
#     h::Real=1.0, bc::SVector{3,<:Real}=SVector(0.0,0.0,0.0)) where {T}

#     ex,ey,ez = m.exp
#     return m.val * face_integral(node_coords, ex, ey, ez, bc) / h^(ex+ey+ez)
# end


# ########## 3D Volume Integration via Divergence Theorem ##########

# """
#     volume_integral(m::Monomial{T,3}, points::AbstractVector{<:SVector{3}},
#                     faces::Vector{<:AbstractVector{Int}};
#                     h::Real=1.0, bc::SVector{3,<:Real}=SVector(0.0,0.0,0.0)) where T

# Compute ∫_V m(((x-bc)/h)) dV over a closed polyhedron given by `points` and `faces`.
# Uses F = ( (x-bcx)^(px+1) (y-bcy)^py (z-bcz)^pz / (px+1), 0, 0 ) / h^(px+py+pz) and sums face contributions.
# """
# function volume_integral(m::Monomial{T,3}, points::AbstractVector{V}, faces::Vector{<:AbstractVector{Int}};
#     h::Real=1.0, bc::SVector{3,<:Real}=SVector(0.0,0.0,0.0), method::Symbol=:faces, tets=nothing) where {T,V<:SVector{3,<:Real}}

#     if method != :faces && method != :tets
#         error("Unknown method=$(method). Use :faces or :tets")
#     end

#     if method == :tets
#         return volume_integral_tets(m, points; h=h, bc=bc, tets=tets)
#     end

#     px,py,pz = m.exp
#     # Interior centroid
#     c = reduce(+, points) / T(length(points))

#     # Prepare monomial exponents for surface integral (x exponent +1)
#     ex_face = px + 1
#     ey_face = py
#     ez_face = pz

#     sum_faces = zero(promote_type(T, eltype(V)))
#     @inbounds for f in faces
#         poly = points[f]
#         # Plane basis and normal
#         o,u,v, _pts2 = project_to_plane(poly)
#         n = cross(u, v) # u,v orthonormal => |n|=1
#         # Ensure outward normal
#         if dot(n, c - poly[1]) > 0
#             n = -n
#         end
#         # Face integral of (x-bc)^(px+1)*(y-bc)^py*(z-bc)^pz
#         If = face_integral(poly, ex_face, ey_face, ez_face, bc)
    
#         sum_faces += n[1] * If
#     end

#     return m.val * (sum_faces / (h^(px+py+pz) * (px+1)))
# end


# ########## 3D Volume Integration via Tetrahedra (Duffy GL product) ##########

# # Optional import of tetrahedralization utility
# try
#     using Main.Triangulation3D: tetrahedralize_points_convex
# catch
#     try
#         import ..Triangulation3D: tetrahedralize_points_convex
#     catch
#         # If not available due to include order, the caller can provide tets explicitly
#     end
# end


# # ========= New: Volume integration via face triangles (pretriangulated) =========

# # Access precomputed face triangulations and topology API without strict include order
# try
#     using Main.Triangulation3D: FaceTriangulations3D, get_area_triangles
# catch
#     try
#         import ..Triangulation3D: FaceTriangulations3D, get_area_triangles
#     catch
#     end
# end

# try
#     using Main: Topology, get_nodes, get_coords, get_volume_node_ids, get_volume_area_ids, get_area_node_ids
# catch
#     try
#         import ..Topology
#         import ..get_nodes
#         import ..get_coords
#         import ..get_volume_node_ids
#         import ..get_volume_area_ids
#         import ..get_area_node_ids
#     catch
#     end
# end

# @inline function _duffy_tri_orders(total_degree::Int)
#     nu = max(1, ceil(Int, (total_degree + 2) / 2))
#     nv = max(1, ceil(Int, (total_degree + 1) / 2))
#     return nu, nv
# end

# @inline function _triangle_flux_x(px::Int, py::Int, pz::Int,
#     p1::StaticVector{3,TS}, p2::StaticVector{3,TS}, p3::StaticVector{3,TS},
#     interior_ref::StaticVector{3,TS}; h::Real=1.0, bc::StaticVector{3,<:Real}=SVector(0.0,0.0,0.0)) where {TS<:Real}

#     total_degree = px + py + pz
#     nu, nv = _duffy_tri_orders(total_degree + 1) # +1 due to px+1 in integrand
#     up, uw = get_gauss_legendre01(nu)
#     vp, vw = get_gauss_legendre01(nv)

#     e1 = p2 - p1
#     e2 = p3 - p1
#     n_area_vec = cross(e1, e2)
#     # Orient outward relative to interior reference
#     if dot(n_area_vec, interior_ref - p1) > 0
#         n_area_vec = -n_area_vec
#     end

#     acc = zero(promote_type(TS, Float64))
#     bx, by, bz = bc
#     n_x = n_area_vec[1]
#     @inbounds for iu in eachindex(up, uw)
#         u = up[iu]
#         wu = uw[iu]
#         one_minus_u = 1 - u
#         @inbounds for iv in eachindex(vp, vw)
#             v = vp[iv]
#             wv = vw[iv]
#             r = u
#             s = v * one_minus_u
#             x = p1 + r*e1 + s*e2
#             # Unscaled monomial for surface integrand (x-bc)^(px+1)*(y-bc)^py*(z-bc)^pz
#             fx = (x[1] - bx)^(px+1) * (x[2] - by)^py * (x[3] - bz)^pz
#             acc += fx * n_x * one_minus_u * wu * wv
#         end
#     end
#     return acc
# end

# """
#     volume_integral_face_tris(m::Monomial{T,3}, topo::Topology{3}, volume_id::Int,
#                                ft::FaceTriangulations3D; h=1.0, bc=SVector(0,0,0)) where T

# Compute ∫_V m(((x-bc)/h)) dV by summing flux over pre-triangulated faces.
# Avoids per-face 3D→2D polynomial projection; integrates directly over 3D triangles
# using Duffy-mapped Gauss–Legendre product.
# """
# function volume_integral_face_tris(m::Monomial{T,3}, topo::Topology{3}, volume_id::Int,
#     ft; h::Real=1.0, bc::SVector{3,<:Real}=SVector(0.0,0.0,0.0)) where {T}

#     px, py, pz = m.exp
#     total_degree = px + py + pz
#     # Interior reference point (simple centroid of volume nodes)
#     vol_gids = get_volume_node_ids(topo, volume_id)
#     vol_points = @views get_nodes(topo)[vol_gids]
#     interior_ref = reduce(+, vol_points) / length(vol_points)

#     sum_flux_x = zero(promote_type(T, Float64))
#     for area_id in get_volume_area_ids(topo, volume_id)
#         gids = get_area_node_ids(topo, area_id)
#         # Access precomputed triangles without relying on imported function names
#         # tris = getproperty(ft, :area_tris)[area_id]
#         tris = ft.area_tris[area_id]
#         isempty(tris) && continue
#         pts = @views get_nodes(topo)[gids]
#         @inbounds for tri in tris
#             i, j, k = tri
#             p1 = pts[i]; p2 = pts[j]; p3 = pts[k]
#             sum_flux_x += _triangle_flux_x(px, py, pz, p1, p2, p3, interior_ref; h=h, bc=bc)
#         end
#     end

#     return m.val * (sum_flux_x / (h^total_degree * (px + 1)))
# end

# # Convenience wrapper with the same name as legacy API
# volume_integral(m::Monomial{T,3}, topo::Topology{3}, volume_id::Int,
#     ft; h::Real=1.0, bc::SVector{3,<:Real}=SVector(0.0,0.0,0.0)) where {T} =
#     volume_integral_face_tris(m, topo, volume_id, ft; h=h, bc=bc)

# @inline function get_gauss_legendre01(n::Int)
#     pts, wts = get_gauss_legendre(n)
#     # Map [-1,1] -> [0,1]
#     pts01 = 0.5 .* (pts .+ 1)
#     wts01 = 0.5 .* wts
#     return pts01, wts01
# end

# @inline function _duffy_tet_orders(total_degree::Int)
#     # Conservative exactness bounds for product GL on [0,1]^3 with Duffy mapping
#     # deg_u ≤ D + 2, deg_v ≤ D + 1, deg_w ≤ D
#     # Use exactness ≈ 2n-1 for Gauss-Legendre; conservative for transformed polynomials
#     nu = max(1, ceil(Int, (total_degree + 3) / 2))
#     nv = max(1, ceil(Int, (total_degree + 2) / 2))
#     nw = max(1, ceil(Int, (total_degree + 1) / 2))
#     return nu, nv, nw
# end

# function tet_volume_integral(m::Monomial{T,3}, p1::SVector{3,TS}, p2::SVector{3,TS}, p3::SVector{3,TS}, p4::SVector{3,TS};
#     h::Real=1.0, bc::SVector{3,<:Real}=SVector(0.0,0.0,0.0)) where {T,TS<:Real}

#     D = m.exp[1] + m.exp[2] + m.exp[3]
#     nu, nv, nw = _duffy_tet_orders(D)
#     up, uw = get_gauss_legendre01(nu)
#     vp, vw = get_gauss_legendre01(nv)
#     wp, ww = get_gauss_legendre01(nw)

#     # Affine map from reference simplex (r,s,t) to physical
#     # r = u, s = v*(1-u), t = w*(1-u)*(1-v)
#     e1 = p2 - p1
#     e2 = p3 - p1
#     e3 = p4 - p1
#     detJ = abs(dot(cross(e1, e2), e3))

#     integral_val = zero(promote_type(T, TS))
#     exx, eyy, ezz = m.exp
#     mval = m.val
#     @inbounds for iu in eachindex(up, uw)
#         u = up[iu]
#         wu = uw[iu]
#         one_minus_u = 1 - u
#         @inbounds for iv in eachindex(vp, vw)
#             v = vp[iv]
#             wv = vw[iv]
#             one_minus_v = 1 - v
#             @inbounds for iw in eachindex(wp, ww)
#                 w = wp[iw]
#                 ww_ = ww[iw]
#                 r = u
#                 s = v * one_minus_u
#                 t = w * one_minus_u * one_minus_v
#                 # Physical point components
#                 x1 = p1[1] + r*e1[1] + s*e2[1] + t*e3[1]
#                 x2 = p1[2] + r*e1[2] + s*e2[2] + t*e3[2]
#                 x3 = p1[3] + r*e1[3] + s*e2[3] + t*e3[3]
#                 # Monomial value at scaled location (avoid SVector allocations)
#                 sx = (x1 - bc[1]) / h
#                 sy = (x2 - bc[2]) / h
#                 sz = (x3 - bc[3]) / h
#                 fx = mval * (sx^exx) * (sy^eyy) * (sz^ezz)
#                 # Jacobian of Duffy mapping (to simplex) and affine map to physical
#                 jac_duffy = one_minus_u^2 * one_minus_v
#                 wprod = wu * wv * ww_
#                 integral_val += fx * jac_duffy * wprod
#             end
#         end
#     end

#     return detJ * integral_val
# end

# function volume_integral_tets(m::Monomial{T,3}, points::AbstractVector{V};
#     h::Real=1.0, bc::SVector{3,<:Real}=SVector(0.0,0.0,0.0), tets=nothing) where {T,V<:SVector{3,<:Real}}

#     local_tets = tets
#     if local_tets === nothing
#         if @isdefined tetrahedralize_points_convex
#             local_tets = tetrahedralize_points_convex(points)
#         else
#             error("No tets provided and tetrahedralize_points_convex not available. Pass tets=... explicitly.")
#         end
#     end

#     acc = zero(promote_type(T, eltype(V)))
#     @inbounds for tet in local_tets
#         a,b,c,d = tet
#         acc += tet_volume_integral(m, points[a], points[b], points[c], points[d]; h=h, bc=bc)
#     end
#     return acc
# end


# ########## Batched monomial basis integrals (up to order O) ##########

# """
#     volume_monomial_moments_tets(points, O; h=1.0, bc=SVector(0,0,0), tets=nothing)

# Compute integrals of all monomials up to total order `O` over a tetrahedral subdivision
# of the polyhedron defined by `points` (and `tets`). Returns a vector `vals` aligned
# with `get_base(BaseInfo{3,O,1}()).base`.

# Uses one quadrature grid sized for the highest degree and accumulates all moments in
# one pass per tet for speed.
# """
# function volume_monomial_moments_tets(points::AbstractVector{V}, O::Int;
#     h::Real=1.0, bc::SVector{3,<:Real}=SVector(0.0,0.0,0.0), tets=nothing) where {V<:SVector{3,<:Real}}

#     # Build/capture basis (U=1)
#     base = get_base(BaseInfo{3,O,1}())
#     L = length(base)
#     Tacc = promote_type(eltype(V), Float64)
#     vals = zeros(Tacc, L)

#     local_tets = tets
#     if local_tets === nothing
#         if @isdefined tetrahedralize_points_convex
#             local_tets = tetrahedralize_points_convex(points)
#         else
#             error("No tets provided and tetrahedralize_points_convex not available. Pass tets=... explicitly.")
#         end
#     end

#     # Quadrature orders for the highest total degree O
#     nu, nv, nw = _duffy_tet_orders(O)
#     up, uw = get_gauss_legendre01(nu)
#     vp, vw = get_gauss_legendre01(nv)
#     wp, ww = get_gauss_legendre01(nw)

#     # Scratch power arrays up to O for each coordinate
#     px = Vector{Tacc}(undef, O+1)
#     py = Vector{Tacc}(undef, O+1)
#     pz = Vector{Tacc}(undef, O+1)

#     @inbounds for tet in local_tets
#         a,b,c,d = tet
#         p1 = points[a]; p2 = points[b]; p3 = points[c]; p4 = points[d]
#         e1 = p2 - p1; e2 = p3 - p1; e3 = p4 - p1
#         detJ = abs(dot(cross(e1, e2), e3))

#         @inbounds for iu in eachindex(up, uw)
#             u = up[iu]; wu = uw[iu]
#             one_minus_u = 1 - u
#             @inbounds for iv in eachindex(vp, vw)
#                 v = vp[iv]; wv = vw[iv]
#                 one_minus_v = 1 - v
#                 @inbounds for iw in eachindex(wp, ww)
#                     w = wp[iw]; ww_ = ww[iw]
#                     r = u
#                     s = v * one_minus_u
#                     t = w * one_minus_u * one_minus_v
#                     x = p1 + r*e1 + s*e2 + t*e3
#                     # scaled coords
#                     sx = (x[1] - bc[1]) / h
#                     sy = (x[2] - bc[2]) / h
#                     sz = (x[3] - bc[3]) / h
#                     # build powers 0..O
#                     px[1] = one(Tacc); py[1] = one(Tacc); pz[1] = one(Tacc)
#                     @inbounds for i in 2:O+1
#                         px[i] = px[i-1] * sx
#                         py[i] = py[i-1] * sy
#                         pz[i] = pz[i-1] * sz
#                     end
#                     # combined weight
#                     w_duffy = one_minus_u^2 * one_minus_v
#                     w_all = detJ * w_duffy * wu * wv * ww_
#                     # accumulate all basis moments
#                     @inbounds for k in 1:L
#                         exp = base.base[k].exp
#                         vals[k] += w_all * px[exp[1]+1] * py[exp[2]+1] * pz[exp[3]+1]
#                     end
#                 end
#             end
#         end
#     end

#     return vals
# end

# """
#     volume_monomial_moments(points, faces, O; method=:tets, h=1.0, bc=SVector(0,0,0), tets=nothing)

# Wrapper that computes moments for all monomials up to `O` either via tets (fast batched)
# or faces (looping over basis; exact but slower).
# """
# function volume_monomial_moments(points::AbstractVector{V}, faces::Vector{<:AbstractVector{Int}}, O::Int;
#     method::Symbol=:tets, h::Real=1.0, bc::SVector{3,<:Real}=SVector(0.0,0.0,0.0), tets=nothing) where {V<:SVector{3,<:Real}}

#     if method == :tets
#         return volume_monomial_moments_tets(points, O; h=h, bc=bc, tets=tets)
#     elseif method == :faces
#         base = get_base(BaseInfo{3,O,1}())
#         L = length(base)
#         vals = zeros(promote_type(eltype(V), Float64), L)
#         @inbounds for k in 1:L
#             exp = base.base[k].exp
#             m = Monomial(1.0, exp)
#             vals[k] = volume_integral(m, points, faces; h=h, bc=bc, method=:faces)
#         end
#         return vals
#     else
#         error("Unknown method=$(method). Use :tets or :faces")
#     end
# end