using StaticArrays, LinearAlgebra, Random
using WriteVTK
using OrderedCollections, Bumper
using Statistics
using SmallCollections, Chairmarks
using JET
using Ju3VEM

################################################################################
# 1) Triangulate a planar polygon in 3D (arbitrary orientation)
#    Export as 2D topology using projected 2D points
################################################################################
let
    # Build a simple concave polygon in 2D
    poly2 = [
        SA[0.0, 0.0],
        SA[1.0, 0.0],
        SA[1.2, 0.5],
        SA[1.0, 1.0],
        SA[0.0, 1.0],
        SA[0.5, 0.5],
    ] |> reverse

    # Random oriented plane via orthonormal basis
    rng = MersenneTwister(1234)
    n = normalize(@SVector randn(rng, 3))
    tmp = abs(n[1]) < 0.9 ? SA[1.0,0.0,0.0] : SA[0.0,1.0,0.0]
    u = normalize(cross(n, tmp))
    v = cross(n, u)
    o = SA[0.3, -0.2, 0.4]

    poly3 = [o + p[1]*u + p[2]*v for p in poly2]

    # Triangulate in 3D
    tris = triangulate_planar_polygon3D(poly3)

    b_res = @b triangulate_planar_polygon3D($poly3)
    display(b_res)

    # Check area sum and positivity
    area_true = polygon_area3D(poly3)
    area_sum = 0.0
    for (i,j,k) in tris
        a,b,c = poly3[i], poly3[j], poly3[k]
        tri_area = 0.5 * norm(cross(b - a, c - a))
        @assert tri_area > 0
        area_sum += tri_area
    end
    @assert isapprox(area_sum, area_true; rtol=1e-10, atol=1e-12)

    # Export as 2D topology using projected 2D points
    o,u,v, pts2 = project_to_plane(poly3)
    topo2 = Topology{2}()
    node_ids2 = Int[]
    append!(node_ids2, add_node!.(pts2, Ref(topo2)))
    for (i,j,k) in tris
        add_area!(SVector(node_ids2[i], node_ids2[j], node_ids2[k]), topo2)
    end
    geometry_to_vtk(topo2, "vtk/triangulated_plane_2d")
end

################################################################################
# 2) Tetrahedralize convex polyhedra (cube and pyramid)
################################################################################
let
    # Cube of side length 1 at origin
    cube_pts = [
        SA[0.0,0.0,0.0], SA[1.0,0.0,0.0], SA[1.0,1.0,0.0], SA[0.0,1.0,0.0],
        SA[0.0,0.0,1.0], SA[1.0,0.0,1.0], SA[1.0,1.0,1.0], SA[0.0,1.0,1.0]
    ]
    topoC = Topology{3}()
    cube_ids = Int[]
    append!(cube_ids, add_node!.(cube_pts, Ref(topoC)))
    faces_cube = [
        [cube_ids[1], cube_ids[2], cube_ids[3], cube_ids[4]], # bottom z=0
        [cube_ids[5], cube_ids[6], cube_ids[7], cube_ids[8]], # top z=1
        [cube_ids[1], cube_ids[2], cube_ids[6], cube_ids[5]], # x edge front
        [cube_ids[2], cube_ids[3], cube_ids[7], cube_ids[6]],
        [cube_ids[3], cube_ids[4], cube_ids[8], cube_ids[7]],
        [cube_ids[4], cube_ids[1], cube_ids[5], cube_ids[8]],
    ]
    add_volume!(faces_cube, topoC)
    vol_id = 1

    tets_local, l2g = tetrahedralize_volume(topoC, vol_id)
    b_res = @b tetrahedralize_volume($topoC, $vol_id)
    display(b_res)
    pts = get_coords.(get_nodes(topoC)[l2g])

    # Check total volume and positivity
    vol_true = 1.0
    vol_sum = 0.0
    for (a,b,c,d) in tets_local
        pa,pb,pc,pd = pts[a], pts[b], pts[c], pts[d]
        v = abs(dot(cross(pb - pa, pc - pa), pd - pa) / 6)
        @assert v > 0
        vol_sum += v
    end
    
    @assert isapprox(vol_sum, vol_true; rtol=1e-6, atol=1e-6)

    # Build a tet-only topology for vtk
    tet_topo = build_tet_topology_from_volume(topoC, vol_id; tets_local)
    

    geometry_to_vtk(tet_topo, "vtk/cube_tets")

    # Pyramid test
    topoP = Topology{3}()
    base = [SA[0.0, 0.0, 0.0], SA[1.0, 0.0, 0.0], SA[1.0, 1.0, 0.0], SA[0.0, 1.0, 0.0]]
    apex = SA[0.4, 0.7, 1.2]
    idsP = Int[]
    append!(idsP, add_node!.(base, Ref(topoP)))
    push!(idsP, add_node!(apex, topoP))
    base_ids = idsP[1:4]; apex_id = idsP[5]
    facesP = [
        base_ids,
        [base_ids[1], base_ids[2], apex_id],
        [base_ids[2], base_ids[3], apex_id],
        [base_ids[3], base_ids[4], apex_id],
        [base_ids[4], base_ids[1], apex_id],
    ]
    add_volume!(facesP, topoP)
    vol_idP = 1

    tets_localP, l2gP = tetrahedralize_volume(topoP, vol_idP)

    ptsP = get_coords.(get_nodes(topoP)[l2gP])

    vol_trueP = volume_of_topo_volume(topoP, vol_idP)
    vol_sumP = 0.0
    for (a,b,c,d) in tets_localP
        pa,pb,pc,pd = ptsP[a], ptsP[b], ptsP[c], ptsP[d]
        v = abs(dot(cross(pb - pa, pc - pa), pd - pa) / 6)
        @assert v > 0
        vol_sumP += v
    end

    @assert isapprox(vol_sumP, vol_trueP; rtol=1e-10, atol=1e-12)

    # tet_topoP = Triangulation3D.build_tet_topology_from_volume(topoP, vol_idP; tets_local=tets_localP)
    # geometry_to_vtk(tet_topoP, "pyramid_tets")
end

################################################################################
# 3) 2D star-like shape (non-self-intersecting), triangulation and export (2D)
################################################################################
let
    # Build a 10-vertex star-shaped decagon by alternating radii in angular order
    R = 1.0
    r = 0.4
    n = 10
    star2d = SVector{2,Float64}[]
    for k in 0:n-1
        θ = 2π * (k / n)
        ρ = isodd(k+1) ? R : r
        push!(star2d, @SVector [ρ*cos(θ), ρ*sin(θ)])
    end

    # Triangulate using 2D ear clipping
    tris = triangulate_polygon(star2d)

    # Polygon area via shoelace
    area_true = 0.0
    for i in 1:length(star2d)
        j = i == length(star2d) ? 1 : i + 1
        area_true += star2d[i][1]*star2d[j][2] - star2d[j][1]*star2d[i][2]
    end
    area_true = abs(area_true) / 2

    # Check area sum and positivity
    area_sum = 0.0
    for (i,j,k) in tris
        a,b,c = star2d[i], star2d[j], star2d[k]
        tri_area = 0.5 * abs((b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1]))
        @assert tri_area > 0
        area_sum += tri_area
    end
    @assert isapprox(area_sum, area_true; rtol=1e-12, atol=1e-12)

    # Export as 2D topology
    topo_star = Topology{2}()
    node_ids = Int[]
    append!(node_ids, add_node!.(star2d, Ref(topo_star)))
    for (i,j,k) in tris
        add_area!(SVector(node_ids[i], node_ids[j], node_ids[k]), topo_star)
    end
    geometry_to_vtk(topo_star, "vtk/star2d_tri")
end

################################################################################
# 4) 3D true polyhedron with polygonal faces: hexagonal prism
################################################################################
let
    m = 6
    R = 1.0
    h = 1.3
    bottom = [@SVector [R*cos(2π*k/m), R*sin(2π*k/m), 0.0] for k in 0:m-1]
    top    = [@SVector [R*cos(2π*k/m), R*sin(2π*k/m), h]   for k in 0:m-1]

    topo = Topology{3}()
    ids_bottom = Int[]
    ids_top = Int[]
    append!(ids_bottom, add_node!.(bottom, Ref(topo)))
    append!(ids_top,    add_node!.(top,    Ref(topo)))

    # Faces: bottom hex, top hex, and 6 quads for sides
    faces = Vector{Vector{Int}}()
    push!(faces, copy(ids_bottom))
    push!(faces, copy(ids_top))
    for i in 1:m
        ip1 = i == m ? 1 : i + 1
        push!(faces, [ids_bottom[i], ids_bottom[ip1], ids_top[ip1], ids_top[i]])
    end

    add_volume!(faces, topo)
    vol_id = 1

    # Tetrahedralize and check volume & positivity
    tets_local, l2g = tetrahedralize_volume(topo, vol_id)
    pts = get_coords.(get_nodes(topo)[l2g])

    # Compute volume by summing oriented volumes, flip if needed per tet
    vol_true = volume_of_topo_volume(topo, vol_id)
    vol_sum = 0.0
    for (a,b,c,d) in tets_local
        pa,pb,pc,pd = pts[a], pts[b], pts[c], pts[d]
        v = dot(cross(pb - pa, pc - pa), pd - pa) / 6
        if v < 0
            # swap orientation
            v = -v
        end
        @assert v > 0
        vol_sum += v
    end
    @assert isapprox(vol_sum, vol_true; rtol=1e-8, atol=1e-10)

    tet_topo = build_tet_topology_from_volume(topo, vol_id; tets_local)
    geometry_to_vtk(tet_topo, "vtk/hex_prism_tets")
end

# ################################################################################
# # 5) 3D polyhedron benchmark: extruded concave arrow-like 2D structure
# ################################################################################
# let
#     # Define a concave arrow polygon (with side overhangs), intended CCW order
#     arrow2d = [
#         @SVector[-1.0, -0.3],
#         @SVector[ 0.3, -0.3],
#         @SVector[ 0.3, -0.6],
#         @SVector[ 1.2,  0.0],
#         @SVector[ 0.3,  0.6],
#         @SVector[ 0.3,  0.3],
#         @SVector[-1.0,  0.3]
#     ]

#     # Ensure CCW orientation for ear clipping
#     s = 0.0
#     for i in 1:length(arrow2d)
#         j = i == length(arrow2d) ? 1 : i + 1
#         s += arrow2d[i][1]*arrow2d[j][2] - arrow2d[j][1]*arrow2d[i][2]
#     end
#     if s < 0
#         arrow2d = reverse(arrow2d)
#     end

#     # Extrude in z to build a non-convex prism
#     h = 1.1
#     bottom = [@SVector [p[1], p[2], 0.0] for p in arrow2d]
#     top    = [@SVector [p[1], p[2], h]   for p in arrow2d]

#     topo = Topology{3}()
#     ids_bottom = Int[]
#     ids_top = Int[]
#     append!(ids_bottom, add_node!.(bottom, Ref(topo)))
#     append!(ids_top,    add_node!.(top,    Ref(topo)))

#     # Build faces: bottom polygon, top polygon, and side quads
#     faces = Vector{Vector{Int}}()
#     push!(faces, copy(ids_bottom))
#     push!(faces, copy(ids_top))
#     m = length(arrow2d)
#     for i in 1:m
#         ip1 = i == m ? 1 : i + 1
#         push!(faces, [ids_bottom[i], ids_bottom[ip1], ids_top[ip1], ids_top[i]])
#     end

#     add_volume!(faces, topo)
#     vol_id = 1

#     # Triangulate the base polygon (handles concavity)
#     tris2d = EarClippingTriangulation2D.triangulate_polygon(arrow2d; find_optimal_ear=true)

#     # Build tetrahedra by decomposing each triangular prism into 3 tets
#     l2g = get_volume_node_ids(topo, vol_id)
#     g2l = Dict(gid => i for (i, gid) in enumerate(l2g))
#     tets_local = NTuple{4,Int}[]
#     for (i,j,k) in tris2d
#         bi, bj, bk = ids_bottom[i], ids_bottom[j], ids_bottom[k]
#         ti, tj, tk = ids_top[i],    ids_top[j],    ids_top[k]
#         # 3-tet decomposition with common apex tk on the top triangle
#         push!(tets_local, (g2l[bi], g2l[bj], g2l[bk], g2l[tk]))
#         push!(tets_local, (g2l[bi], g2l[bj], g2l[tj], g2l[tk]))
#         push!(tets_local, (g2l[bi], g2l[ti], g2l[tj], g2l[tk]))
#     end

#     # Benchmark the base triangulation
#     b_res = @b EarClippingTriangulation2D.triangulate_polygon($arrow2d)
#     display(b_res)

#     # Validate volume & positivity
#     pts = get_coords.(get_nodes(topo)[l2g])
#     vol_true = Triangulation3D.volume_of_topo_volume(topo, vol_id)
#     vol_sum = 0.0
#     for (a,b,c,d) in tets_local
#         pa,pb,pc,pd = pts[a], pts[b], pts[c], pts[d]
#         v = abs(dot(cross(pb - pa, pc - pa), pd - pa) / 6)
#         @assert v > 0
#         vol_sum += v
#     end
#     @show vol_sum, vol_true
#     # @assert isapprox(vol_sum, vol_true; rtol=1e-8, atol=1e-12)

#     tet_topo = Triangulation3D.build_tet_topology_from_volume(topo, vol_id; tets_local)
#     geometry_to_vtk(tet_topo, "arrow_tets")
# end

nothing


