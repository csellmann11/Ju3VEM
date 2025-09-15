"""
    FaceTriangulations3D

Precompute and store a triangulation for every face (area) in a `Topology{3}`.
Stores triangles as local indices into each face's node list.
"""
struct FaceTriangulations3D
    topo::Topology{3}
    area_tris::Vector{Vector{NTuple{3,Int}}}
end

function FaceTriangulations3D(topo::Topology{3}; find_optimal_ear::Bool=false)
    nareas = length(get_areas(topo))
    tris_per_area = Vector{Vector{NTuple{3,Int}}}(undef, nareas)
    for area_id in 1:nareas
        gids = get_area_node_ids(topo, area_id)
        if isempty(gids)
            tris_per_area[area_id] = NTuple{3,Int}[]
            continue
        end
        local_tris = triangulate_area_local_ids(topo, area_id; find_optimal_ear=find_optimal_ear)
        tris = Vector{NTuple{3,Int}}(undef, length(local_tris))
        for (i,t) in enumerate(local_tris)
            tris[i] = (t[1], t[2], t[3])
        end
        tris_per_area[area_id] = tris
    end
    return FaceTriangulations3D(topo, tris_per_area)
end


struct FaceIntegrals{D,L} 

    integrals::Vector{MVector{L,Float64}}
    function FaceIntegrals{D}(integrals::Vector{MVector{L,Float64}}) where {D,L}
        new{D,L}(integrals)
    end
end


function get_face_integral(exp::StaticVector{D,Int},face_id::Int,fi::FaceIntegrals{D,L}) where {D,L}
    idx = get_exp_to_idx_dict(exp)
    return fi.integrals[face_id][idx]
end

function get_face_integral(m::Monomial,face_id::Int,fi::FaceIntegrals{D,L}) where {D,L}
    return get_face_integral(m.exp,face_id,fi)
end


function precompute_face_integrals(topo::Topology{D}, 
    ft::FaceTriangulations3D,bi::BaseInfo{D,O}) where {D,O}
    nareas = length(get_areas(topo))
    base      = get_base(bi)
    integrals = Vector{MVector{length(base),Float64}}(undef, nareas)
    for area in RootIterator{D,3}(topo)
        for (i,m) in enumerate(base.base)
            integrals[area.id][i] = integrate_polynomial_over_face(m, area.id, topo, ft)
        end
    end
    return FaceIntegrals{D}(integrals)
end




@inline get_area_triangles(ft::FaceTriangulations3D, area_id::Int) = ft.area_tris[area_id]

function integrate_polynomial_over_face(m::Monomial, face_id::Int, topo::Topology{3}, ft)
    area_node_ids = get_area_node_ids(topo, face_id)
    tris = get_area_triangles(ft, face_id)

    n = get_plane_parameters(@views get_nodes(topo)[area_node_ids])[4]
    order = min(1, sum(m.exp))
    quad_rule = _TriangleQuadRuleLookup[order]

    int = 0.0
    for tri in tris
        i,j,k = tri
        p1 = get_nodes(topo)[area_node_ids[i]]
        p2 = get_nodes(topo)[area_node_ids[j]]
        p3 = get_nodes(topo)[area_node_ids[k]]
        nx = n[1] 

        for (w,p) in zip(quad_rule.weights, quad_rule.points)
            ξ1,ξ2 = p
            ξ3    = 1-ξ1-ξ2 
            x     = p1*ξ1 + p2*ξ2 + p3*ξ3
            int  += w*m(x)*nx
        end
    end
    return int, n
end

function integrate_polynomial_over_volume(m::Monomial, tet_id::Int, topo::Topology{3}, ft)

    exp = SVector(m.exp[1]+1, m.exp[2], m.exp[3])
    quad_rule = _TriangleQuadRuleLookup[sum(exp)]

    mface       = Monomial(m.val, exp)
    v_node_ids  = get_volume_node_ids(topo, tet_id)

    #! make this more robust
    element_center = mean(get_nodes(topo)[vi] for vi in v_node_ids)

    int = 0.0
    for area_id in get_volume_area_ids(topo, tet_id)
        _area_int,n = integrate_polynomial_over_face(mface, area_id, topo, ft) 
        face_node   = get_nodes(topo)[get_area_node_ids(topo, area_id)[1]]
        λ           = check_normal_sign(n, face_node, element_center)
        int        += λ*_area_int
    end 
    # return int/(exp[1])
    int/(exp[1])
end
