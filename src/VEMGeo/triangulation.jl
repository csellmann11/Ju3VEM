using StaticArrays
using Random
using LinearAlgebra
using FixedSizeArrays

"""
    triangle_quality_measure(p1::AbstractVector{V}, p2::AbstractVector{V}, p3::AbstractVector{V}) where {V<:AbstractVector}

Calculates a quality measure for a triangle defined by three nodes.
The measure is based on the ratio of the triangle's area to the sum of the squared side lengths.
A higher value indicates a better-shaped triangle (closer to equilateral).

This measure is normalized such that an equilateral triangle has a value of approximately 0.577.

# Arguments
- `p1`, `p2`, `p3`: Three `AbstractVector`s representing the coordinates of the triangle's nodes.
                    Each `AbstractVector` should have the same dimension (2D or 3D).

# Returns
- A `Float64` value representing the quality measure.
"""
function triangle_quality_measure(p1::AbstractVector{T}, p2::AbstractVector{T}, p3::AbstractVector{T}) where {T<:Real}
    # Calculate side vectors
    v21 = p2 - p1
    v31 = p3 - p1
    v32 = p3 - p2

    # Calculate squared side lengths
    a_sq = dot(v32, v32) # |p3 - p2|^2
    b_sq = dot(v31, v31) # |p3 - p1|^2
    c_sq = dot(v21, v21) # |p2 - p1|^2

    # Calculate area (using cross product for 3D, or 2D equivalent)
    # For 2D points (represented as 2-element vectors), a 2D cross product equivalent is used.
    # For 3D points, the magnitude of the 3D cross product is 2 * Area.
    area = if length(p1) == 2
        0.5 * abs(p1[1] * (p2[2] - p3[2]) + p2[1] * (p3[2] - p1[2]) + p3[1] * (p1[2] - p2[2]))
    elseif length(p1) == 3
        0.5 * norm(cross(v21, v31))
    else
        error("Unsupported dimension for triangle nodes. Only 2D or 3D are supported.")
    end

    # Calculate quality measure
    # The factor of 4*sqrt(3) normalizes the measure for an equilateral triangle.
    # For an equilateral triangle with side 's', Area = (sqrt(3)/4)*s^2, SumSqSides = 3*s^2
    # Quality = Area / SumSqSides = (sqrt(3)/4)*s^2 / (3*s^2) = sqrt(3)/12 approx 0.144
    # Multiplying by 4*sqrt(3) for normalization: (sqrt(3)/12) * 4*sqrt(3) = (4*3)/12 = 1.0
    # A common quality measure is 4 * sqrt(3) * Area / (a^2 + b^2 + c^2)
    # This gives 1.0 for equilateral triangles.
    # Or simply Area / (a_sq + b_sq + c_sq) for a relative measure.
    # Let's use the normalized one that gives 1.0 for equilateral.
    if area == 0.0 || (a_sq + b_sq + c_sq) == 0.0
        return 0.0 # Degenerate triangle
    end

    return (4 * sqrt(3) * area) / (a_sq + b_sq + c_sq)
end


"""
    triangulate_polygon(nodes::AbstractVector{Vec{2,T}}) where T

Triangulate a simple polygon using the ear clipping algorithm.

# Arguments
- `nodes`: Vector of 2D points representing the polygon vertices in order

# Returns
- `triangles`: Vector of tuples containing indices of triangle vertices
"""

function triangulate_polygon(nodes::AbstractVector{V};
    find_optimal_ear::Bool = false) where {V<:AbstractVector{<:Real}}
    n = length(nodes)

    D = 2 # assumes statically sized vector
    @assert D == length(V)

    # Need at least 3 vertices for a triangle
    if n < 3
        return Vector{Tuple{Int,Int,Int}}()
    elseif n == 3
        return [(1, 2, 3)]
    end


    # Create a list of vertex indices
    vertex_indices = collect(1:n)
    # triangles = Vector{Tuple{Int,Int,Int}}()
    n_triangles = length(nodes) - 2
    triangles = FixedSizeArray{NTuple{3,Int}}(undef, n_triangles)
    tria_counter = 1

    # Remove collinear vertices first (optional optimization)
    # vertex_indices = filter_collinear(nodes, vertex_indices)
    if V <:Node
        ids = get_id.(nodes)
        rng = MersenneTwister(hash(ids))
        offset = rand(rng, 1:length(nodes))
    else 
        offset = 1
    end
 
    compare_fun = find_optimal_ear ? !isless : isless
    
    while length(vertex_indices) > 3
        n_current = length(vertex_indices)
        ear_found = false

        tria_measure = find_optimal_ear ? -Inf : Inf
        best_triangle = (tria_measure, (-1, -1, -1), -1)
        # Try to find a valid ear
        # for i in 1:n_current
        for _i in offset:(n_current+offset-1)
            i = mod1(_i, n_current)
            prev_i = mod1(i - 1, n_current)
            next_i = mod1(i + 1, n_current)

            # Get actual vertex indices
            prev_idx = vertex_indices[prev_i]
            curr_idx = vertex_indices[i]
            next_idx = vertex_indices[next_i]

            # Get vertex positions
            v_prev = SVector{D}(nodes[prev_idx])
            v_curr = SVector{D}(nodes[curr_idx])
            v_next = SVector{D}(nodes[next_idx])

            # Check if this forms a valid ear
            if is_valid_ear(v_prev, v_curr, v_next, nodes, vertex_indices, prev_i, i, next_i)

                tria_measure = triangle_quality_measure(v_prev, v_curr, v_next)

                if compare_fun(tria_measure,best_triangle[1])
                # if tria_measure > best_triangle[1]
                    best_triangle = (tria_measure, (prev_idx, curr_idx, next_idx), i)
                end
                ear_found = true
                !find_optimal_ear && break
            end
        end

        triangles[tria_counter] = best_triangle[2]
        tria_counter += 1

        if !ear_found throw(ErrorException("No ear found")) end
        deleteat!(vertex_indices, best_triangle[3])
    end

    # Add the final triangle
    if length(vertex_indices) == 3
        triangles[end] = (vertex_indices[1], vertex_indices[2], vertex_indices[3])
    end

    return triangles
end

"""
    is_valid_ear(v_prev, v_curr, v_next, nodes, vertex_indices, prev_i, curr_i, next_i)

Check if a triangle forms a valid ear that can be clipped.
"""
function is_valid_ear(v_prev::V, v_curr::V, v_next::V,
    nodes::AbstractVector{<:AbstractVector{<:Real}},
    vertex_indices::Vector{Int},
    prev_i::Int, curr_i::Int, next_i::Int) where {V<:AbstractVector{<:Real}}

    D = 2 # assumes statically sized vector
    @assert D == length(V)

    # First check if the triangle has a positive area (is CCW)
    area = signed_area_2x(v_prev, v_curr, v_next)
    if area <= sqrt(eps(Float64))
        return false  # Reflex vertex or degenerate triangle
    end

    # Check if any other vertex is inside this triangle
    n = length(vertex_indices)
    for i in 1:n
        if i == prev_i || i == curr_i || i == next_i
            continue
        end

        test_vertex = SVector{D}(nodes[vertex_indices[i]])
        if is_point_in_triangle(test_vertex, v_prev, v_curr, v_next)
            return false
        end
    end

    return true
end




"""
    signed_area_2x(v1, v2, v3)

Calculate twice the signed area of a triangle.
Positive if CCW, negative if CW.
"""
function signed_area_2x(v1::V1, v2::V2, v3::V3) where {
    V1<:AbstractVector{<:Real}, V2<:AbstractVector{<:Real}, V3<:AbstractVector{<:Real}}

    return (v2[1] - v1[1]) * (v3[2] - v1[2]) - (v3[1] - v1[1]) * (v2[2] - v1[2])
end

"""
    is_point_in_triangle(p, v1, v2, v3)

Check if point p is strictly inside triangle (v1, v2, v3).
Uses barycentric coordinates with proper handling of edge cases.
"""
function is_point_in_triangle(p::V1, v1::V2, v2::V3, v3::V4) where {
    V1<:AbstractVector{<:Real}, V2<:AbstractVector{<:Real}, V3<:AbstractVector{<:Real}, V4<:AbstractVector{<:Real}}

    D = length(V1)
    # Calculate barycentric coordinates
    v0 = v3 - v1
    v1_vec = v2 - v1
    v2_vec = p - v1

    dot00 = dot(v0, v0)
    dot01 = dot(v0, v1_vec)
    dot02 = dot(v0, v2_vec)
    dot11 = dot(v1_vec, v1_vec)
    dot12 = dot(v1_vec, v2_vec)

    # Calculate barycentric coordinates
    denom = dot00 * dot11 - dot01 * dot01

    # Check for degenerate triangle
    T = promote_type(eltype(p), eltype(v1), eltype(v2), eltype(v3))
    if abs(denom) < eps(T)
        return false
    end

    inv_denom = 1 / denom
    u = (dot11 * dot02 - dot01 * dot12) * inv_denom
    v = (dot00 * dot12 - dot01 * dot02) * inv_denom

    # Check if point is inside Triangle 
    tol = -sqrt(eps(T))
    return (u ≥ tol) && (v ≥ tol) && (u + v ≤ 1 + tol)
end

# Alternative implementation using cross products
function is_point_in_triangle_cross(p::StaticVector{2,T}, v1::StaticVector{2,T}, v2::StaticVector{2,T}, v3::StaticVector{2,T}) where T
    # Check on which side of each edge the point lies
    d1 = sign_of_point_to_line(p, v1, v2)
    d2 = sign_of_point_to_line(p, v2, v3)
    d3 = sign_of_point_to_line(p, v3, v1)

    # Point is inside if all signs are the same
    has_neg = (d1 <= 0) || (d2 <= 0) || (d3 <= 0)
    has_pos = (d1 >= 0) || (d2 >= 0) || (d3 >= 0)

    return !(has_neg && has_pos)
end

"""
    sign_of_point_to_line(p, a, b)

Determine which side of line (a,b) point p is on.
Returns positive if p is to the left, negative if to the right.
"""
function sign_of_point_to_line(p::V1, a::V2, b::V3) where {
    V1<:AbstractVector{<:Real}, V2<:AbstractVector{<:Real}, V3<:AbstractVector{<:Real}}

    D = 2 # assumes statically sized vector
    @assert D == length(V1)
    @assert D == length(V2)
    @assert D == length(V3)
    return (b[1] - a[1]) * (p[2] - a[2]) - (b[2] - a[2]) * (p[1] - a[1])
end

using Statistics

# Quadrature for triangle integration
include("gauss_quad_tri_table.jl")

const Vec3{T} = SVector{3,T}

function get_plane_parameters(points::AbstractVector{<:StaticVector{3,T}}) where T

    tol = sqrt(eps(T))
    o = points[1]

    e1 = points[2] - o
    # Find point that is not collinear with e1
    n = zero(Vec3{T})
    i3 = 3 
    while i3 ≤ length(points)
        e2 = points[i3] - o
        n = cross(e1, e2)
        norm(n) > tol && break
        i3 += 1
    end
    @assert i3 ≤ length(points) "Points appear to be collinear"
    nscaled = normalize(n)
    u = normalize(e1 - dot(e1, nscaled) * nscaled)
    v = cross(nscaled, u)

    return u,v,nscaled,n,o
end



"""
    project_to_plane(points::AbstractVector{<:SVector{3,T}}) where T

Compute an orthonormal basis for the plane spanned by `points` and project them to 2D.

Returns `(origin, u, v, projected_points)`, where `origin` is on the plane,
`u` and `v` are orthonormal in-plane basis vectors, and `projected_points` is
`Vector{SVector{2,T}}`.
"""
function project_to_plane!(points2::AbstractVector{<:StaticVector{2,T}}, 
    points::AbstractVector{<:StaticVector{3,T}}) where T
    @assert length(points) ≥ 3 "Need at least 3 points"
    @assert length(points2) == length(points)

    u,v,n,_,_ = get_plane_parameters(points)
    p0      = points[1]

    for (i,p) in enumerate(points)
        d = p - p0
        points2[i] = SVector(dot(d, u), dot(d, v))
    end
    return u, v, n
end



function project_to_plane(points::AbstractVector{<:SVector{3,T}}) where T
    points2 = Vector{SVector{2,T}}(undef, length(points))
    u,v,n = project_to_plane!(points2, points)
    return u,v,n,points2
end

"""
    triangulate_planar_polygon3D(points::AbstractVector{<:SVector{3,T}}; find_optimal_ear=true) where T

Triangulate a planar polygon in 3D with arbitrary orientation. Returns a vector
of triangle index tuples referring to the input order.
The algorithm projects all points on the plane in the main direction of the normal vector. 
The projection simply drops the component of the point in the main direction of the normal vector.
"""
function triangulate_planar_polygon3D(points::AbstractVector{<:StaticVector{3,T}}; find_optimal_ear::Bool=false) where T
    n = length(points)

    u,v,normal,_,_ = get_plane_parameters(points)
    zero_dir = argmax(abs.(normal))

    pts2 = FixedSizeVector{SVector{2,T}}(undef, n)
    if zero_dir == 1 
        for i in eachindex(points)
            pts2[i] = SVector(points[i][2], points[i][3])
        end
    elseif zero_dir == 2
        for i in eachindex(points)
            pts2[i] = SVector(points[i][1], points[i][3])
        end 
    elseif zero_dir == 3
        for i in eachindex(points)
            pts2[i] = SVector(points[i][1], points[i][2])
        end
    end

    s = zero(T)
    for i in eachindex(pts2)
        j = i == n ? 1 : i + 1
        s += pts2[i][1]*pts2[j][2] - pts2[j][1]*pts2[i][2]
    end

    pts2o = s<0 ? reverse(pts2) : pts2

    return triangulate_polygon(pts2o; find_optimal_ear)
end

@inline signed_tet_volume(a::StaticVector{3,T}, b::StaticVector{3,T}, c::StaticVector{3,T}, d::StaticVector{3,T}) where T =
    dot(cross(b - a, c - a), d - a) / T(6)
@inline tet_volume(a::StaticVector{3,T}, b::StaticVector{3,T}, c::StaticVector{3,T}, d::StaticVector{3,T}) where T =
    abs(signed_tet_volume(a,b,c,d))

"""
    tetrahedralize_points_convex(points::AbstractVector{<:SVector{3,T}}) where T

Given a set of 3D points that are vertices of a convex polyhedron, return a list
of tetrahedra as tuples of local point indices that triangulate the convex hull
without introducing new points. Uses an incremental beneath-beyond style algorithm.
"""
function tetrahedralize_points_convex(points::AbstractVector{<:StaticVector{3,T}}) where T
    N = length(points)
    @assert N ≥ 4 "Need at least 4 non-coplanar points"

    # Interior reference from seed tetra (remains inside as hull grows)

    # Find 4 non-coplanar points to seed
    i1 = 1
    i2 = 2 # assumes that all points are unique
    # pick i3 not collinear
    i3 = -1
    for i in 3:N
        if i != i1 && i != i2
            n = cross(points[i2] - points[i1], points[i] - points[i1])
            if norm(n) > sqrt(eps(T))
                i3 = i
                break
            end
        end
    end
    @assert i3 != -1 "Points appear to be collinear"
    i3 = i3::Int
    # pick i4 not coplanar
    i4 = -1
    for i in 4:N
        if i != i1 && i != i2 && i != i3
            vol = signed_tet_volume(points[i1], points[i2], points[i3], points[i])
            if abs(vol) > sqrt(eps(T))
                i4 = i
                break
            end
        end
    end
    @assert i4 != -1 "Points appear to be coplanar"

    used = SA[i1,i2,i3,i4]
    # Ensure positive orientation for seed tet
    a,b,c,d = i1,i2,i3,i4
    vol = signed_tet_volume(points[a], points[b], points[c], points[d])
    if vol < 0
        c,d = d,c
    end

    # Boundary faces: vector of oriented triples, outward normals oriented away from interior_ref
    function orient_face_outward(tri::NTuple{3,Int}, interior_ref::SVector{3,T})
        i,j,k = tri
        n = cross(points[j] - points[i], points[k] - points[i])
        if dot(n, interior_ref - points[i]) > 0
            # normal points towards centroid; flip to make it outward
            return (i,k,j)
        else
            return tri
        end
    end

    #TODO: Replace this in the future with the barycenter
    interior_ref = (points[a] + points[b] + points[c] + points[d]) / T(4)
   
    boundary = FixedSizeArray{NTuple{3,Int}}(undef, 4)
    boundary[1] = orient_face_outward((a,b,c), interior_ref)
    boundary[2] = orient_face_outward((a,d,b), interior_ref)
    boundary[3] = orient_face_outward((b,d,c), interior_ref)
    boundary[4] = orient_face_outward((c,d,a), interior_ref)
  

    #! This bound might be to small
    max_tet_number = 2*length(points) - 7
    tet_counter = 1
    tets = FixedSizeArray{NTuple{4,Int}}(undef, max_tet_number)

    new_boundary = NTuple{3,Int}[]
    horizon_oriented = Tuple{Int,Int}[]
    tets[1] = (a,b,c,d)
    tet_counter += 1

    edge_map = Dict{Tuple{Int,Int}, Tuple{Int,Tuple{Int,Int}}}()
    for p in 1:N 
        p in used && continue
        # Determine visible faces of boundary from p
        visible = FixedSizeArray{Bool}(undef, length(boundary))
        for (idx, (i,j,k)) in enumerate(boundary)
            n = cross(points[j] - points[i], points[k] - points[i])
            if dot(n, points[p] - points[i]) > sqrt(eps(T))
                visible[idx] = true
            else 
                visible[idx] = false
            end
        end
        if !any(visible)
            # p lies inside current hull due to numeric noise; skip
            continue
        end
        # Add tets from p to all visible faces
        for (idx, vis) in enumerate(visible)
            if vis
                i,j,k = boundary[idx]
                tet = (i,j,k,p)
                # Ensure positive orientation
                vol = signed_tet_volume(points[i], points[j], points[k], points[p])
                if vol < 0
                    tet = (i,k,j,p)
                end
                tets[tet_counter] = tet
                tet_counter += 1
            end
        end
        # Collect horizon edges = edges on boundary between visible and non-visible regions
        # Map unordered edge => (count_visible, one_oriented_edge)
        # edge_map = Dict{Tuple{Int,Int}, Tuple{Int,Tuple{Int,Int}}}()
        empty!(edge_map)
        for (idx, (i,j,k)) in enumerate(boundary)
            if visible[idx]
                for (a1,b1) in ((i,j),(j,k),(k,i))
                    key = a1 < b1 ? (a1,b1) : (b1,a1)
                    entry = get(edge_map, key, (0, (a1,b1)))
                    edge_map[key] = (entry[1] + 1, entry[2])
                end
            end
        end
        empty!(horizon_oriented)
        for (cnt, oriented) in values(edge_map)
            if cnt == 1
                push!(horizon_oriented, oriented)
            end
        end
        # Build new boundary: remove visible faces, add new faces (p with horizon edges)
        empty!(new_boundary)
        for (idx, tri) in enumerate(boundary)
            if !visible[idx]
                push!(new_boundary, tri)
            end
        end
        for (a1,b1) in horizon_oriented
            tri = (a1,b1,p)
            # Orient outward
            n = cross(points[tri[2]] - points[tri[1]], points[tri[3]] - points[tri[1]])
            if dot(n, interior_ref - points[a1]) > 0
                tri = (a1,tri[3],b1)
            end
            push!(new_boundary, tri)
        end
        # boundary = new_boundary
        boundary = FixedSizeArray(new_boundary)
    end

    return tets[1:tet_counter-1]
end

"""
    tetrahedralize_volume_local_ids(topo::Topology{3}, volume_id::Int)

Return tetrahedra as local indices for the given `volume_id` in `topo`.
"""
function tetrahedralize_volume_local_ids(topo::Topology{3}, volume_id::Int)
    global_node_ids = get_volume_node_ids(topo, volume_id)
    # pts = @views get_coords.(get_nodes(topo)[global_node_ids])
    pts = @view get_nodes(topo)[global_node_ids]
    tets_local = tetrahedralize_points_convex(pts)
    return tets_local
end

"""
    tetrahedralize_volume(topo::Topology{3}, volume_id::Int)

Return `(tets_local, local_to_global)` for the given volume.
"""
function tetrahedralize_volume(topo::Topology{3}, volume_id::Int)
    global_node_ids = get_volume_node_ids(topo, volume_id)
    # pts = @views get_coords.(get_nodes(topo)[global_node_ids])
    pts = @view get_nodes(topo)[global_node_ids]
    tets_local = tetrahedralize_points_convex(pts)
    return tets_local, global_node_ids
end

"""
    build_tet_topology_from_volume(topo::Topology{3}, volume_id::Int; tets_local)

Build a new `Topology{3}` containing tetrahedra for visualization. Reuses the
coordinates of the original `topo` and maps local indices to global ones.
"""
function build_tet_topology_from_volume(topo::Topology{3}, volume_id::Int; tets_local)
    tet_topo = Topology{3}()
    local_to_global = get_volume_node_ids(topo, volume_id)
    # Add nodes in local order
    for gid in local_to_global
        add_node!(get_coords(get_nodes(topo)[gid]), tet_topo)
    end
    # Helper to add a tetrahedron as a polyhedron with 4 triangular faces
    function add_tet!(a::Int,b::Int,c::Int,d::Int)
        faces = (
            SVector(a,b,c),
            SVector(a,d,b),
            SVector(b,d,c),
            SVector(c,d,a),
        )
        face_ids = Int[]
        for f in faces
            push!(face_ids, add_area!(f, tet_topo))
        end
        add_volume!(face_ids, tet_topo)
    end
    for tet in tets_local
        add_tet!(tet...)
    end
    return tet_topo
end

"""
    triangulate_area_local_ids(topo::Topology{3}, area_id::Int; find_optimal_ear=true)

Return triangles as local indices into the area's node list.
"""
function triangulate_area_local_ids(topo::Topology{3}, area_id::Int; find_optimal_ear::Bool=true)
    gids = get_area_node_ids(topo, area_id)
    pts = get_nodes(topo)[gids]
    tris_local = triangulate_planar_polygon3D(pts; find_optimal_ear)
    return tris_local
end

"""
    triangulate_area(topo::Topology{3}, area_id::Int; find_optimal_ear=true)

Return triangles as global node id tuples for the given area.
"""
function triangulate_area(topo::Topology{3}, area_id::Int; find_optimal_ear::Bool=true)
    gids = get_area_node_ids(topo, area_id)
    tris_local = triangulate_area_local_ids(topo, area_id; find_optimal_ear)
    return map(t -> (gids[t[1]], gids[t[2]], gids[t[3]]), tris_local)
end



"""
    polygon_area3D(points::AbstractVector{<:SVector{3,T}}) where T

Signed area (positive) of a simple planar polygon in 3D. Only needed for testing.
"""
function polygon_area3D(points::AbstractVector{<:StaticVector{3,T}}) where T
    _,_,_, pts2 = project_to_plane(points)
    # Shoelace formula
    area2 = zero(T)
    for i in eachindex(pts2)
        j = i == length(pts2) ? 1 : i+1
        area2 += pts2[i][1]*pts2[j][2] - pts2[j][1]*pts2[i][2]
    end
    return abs(area2) / 2
end

"""
    volume_from_faces(points::AbstractVector{<:SVector{3,T}}, faces::Vector{<:AbstractVector{Int}}) where T

Compute volume of a closed polyhedron given vertex coordinates and faces (each a
polygon as a list of vertex indices, outward orientation unknown). Triangulates
each face robustly in its own plane (handles concave faces) and sums oriented
tetra volumes via the divergence theorem; triangle orientations are flipped to
be outward using the polyhedron centroid as reference.
"""
function volume_from_faces(points::AbstractVector{<:StaticVector{3,T}}, faces::Vector{<:AbstractVector{Int}}) where T
    c = reduce(+, points) / T(length(points))
    V = zero(T)
    for f in faces
        poly = points[f]
        # Triangulate the face polygon in 3D (handles concavity)
        tris = triangulate_planar_polygon3D(poly; find_optimal_ear=true)
        # Average face normal from triangles (based on current orientation)
        n_avg = zero(SVector{3,T})
        for (i,j,k) in tris
            n_avg += cross(poly[j] - poly[i], poly[k] - poly[i])
        end
        # If normal points toward centroid, flip triangle orientation
        flip = dot(n_avg, c - poly[1]) > 0
        if flip
            for (i,j,k) in tris
                a = poly[i]; b = poly[k]; d = poly[j]
                V += dot(a, cross(b, d)) / T(6)
            end
        else
            for (i,j,k) in tris
                a = poly[i]; b = poly[j]; d = poly[k]
                V += dot(a, cross(b, d)) / T(6)
            end
        end
    end
    return abs(V)
end

"""
    volume_of_topo_volume(topo::Topology{3}, volume_id::Int)

Compute volume of a `Topology` volume using its faces.
"""
function volume_of_topo_volume(topo::Topology{3}, volume_id::Int)
    node_ids = get_volume_node_ids(topo, volume_id)
    pts = get_nodes(topo)[node_ids]
    face_ids = get_volume_area_ids(topo, volume_id)
    faces = get_area_node_ids(topo)[face_ids]
    return volume_from_faces(pts, map(f -> map(x -> findfirst(==(x), node_ids), f), faces))
end


function check_normal_sign(normal, face_node, element_center)
    if dot(normal, face_node - element_center) > 0
        return 1
    else
        return -1
    end
end

function triangle_area_3d_and_normal(A, B, C)
    # Compute vectors from A to B and A to C
    AB = B - A
    AC = C - A
    
    # Area = 0.5 * magnitude of cross product
    unscaled_normal = cross(AB, AC)
    return 0.5 * norm(unscaled_normal), normalize(unscaled_normal)
end




 

