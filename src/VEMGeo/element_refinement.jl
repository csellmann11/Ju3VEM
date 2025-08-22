using SmallCollections

function _refine!(edge::Edge{D},topo::Topology{D}) where D
    edge_id = edge.id

    if !edge.is_root
        mid_node_id = find_single_intersec(get_edge_node_ids(topo)[edge.childs])
        return mid_node_id
    end
 
    edge.is_root = false
    edge_nodes = get_edge_node_ids(topo,edge_id)

    was_once_refined = !isempty(edge.childs)

    if !was_once_refined
        new_node_coords = mean(get_nodes(topo)[n_id] for n_id in edge_nodes)
        new_node_id = add_node!(new_node_coords,topo)
    else 
        new_node_id = find_single_intersec(get_edge_node_ids(topo)[edge.childs])
        activate!(get_nodes(topo)[new_node_id])
        for child_id in edge.childs # reactivate childs
            activate!(get_edges(topo)[child_id])
        end
        return new_node_id
    end

     
    eid1 = add_edge!(SA[edge_nodes[1],new_node_id],
            topo,edge_id,edge.refinement_level)
    eid2 = add_edge!(SA[new_node_id,edge_nodes[2]],
            topo,edge_id,edge.refinement_level)

    push!(edge.childs,eid1); push!(edge.childs,eid2)

    return new_node_id
end


"""
    _refine!(area::Area{D},topo::Topology{D}) where D

Refines the area by adding a new node in the barycenter of the area. The new elements 
do consist out of 4 geometrical nodes. If the edges had hanging nodes, the new elements
will have those nodes as well. Refinement is only save, if element 
is starshaped w.r.t. the barycenter.

# !!! Warning
- Dont use this function directly, use ```refine!(mesh::Mesh{D},elt::Type{ElT},
    marker::BitVector = BitVector(
        1 for _ in get_geo_vec(mesh.topo,Val(D+1)))) where {D,ElT<:ElType} instead```

    -> This function does also update the neccessary informations for the mesh

    
# Arguments
- area::Area{D}: The area to refine
- topo::Topology{D}: The topology of the mesh

# Returns
- new_node_id::Int: The id of the new node at the barycenter
"""
function _refine!(area::Area{D},topo::Topology{D}) where D
    area_id = area.id

    if !is_root(area)
        mid_node = find_single_intersec(get_area_node_ids(topo),area.childs)
        return mid_node
    end

    area.is_root = false
    area_nodes = get_area_node_ids(topo,area_id)

    for edgeid in get_area_edge_ids(topo,area_id)
        edge = get_edges(topo)[edgeid] 
        _refine!(edge,topo)
    end
    
    was_once_refined = !isempty(area.childs)
    if !was_once_refined #standard case
        mid_node_coords = mean(get_nodes(topo)[n_id] for n_id in area_nodes)
        new_node_id     = add_node!(mid_node_coords,topo)
    else
        new_node_id = find_single_intersec(get_area_node_ids(topo),area.childs)
        activate!(get_nodes(topo)[new_node_id])
        # reactivate childs and inner edges
        for child_id in area.childs
            activate!(get_areas(topo)[child_id])
            for edge_id in get_area_edge_ids(topo,child_id)
                activate!(get_edges(topo)[edge_id])
            end
        end
        return new_node_id
    end



    new_geo_nodes = get_iterative_area_vertex_ids(area,topo,
            area.refinement_level+1)

    # Create quadrilaterals by connecting barycenter to refined edge nodes
    for i in eachindex(area_nodes)
        # For first sub-area: use last node, otherwise use node at 2i-2
        first_node = i == 1 ? new_geo_nodes[end] : new_geo_nodes[2i-2] 
        
        # Create quadrilateral: [first_node, edge_node1, edge_node2, barycenter]
        quadrilateral_nodes = SVector(first_node, new_geo_nodes[2i-1], new_geo_nodes[2i], new_node_id)

        narea_id = add_area!(quadrilateral_nodes, topo, area_id, area.refinement_level)
        push!(area.childs, narea_id)
    end
    return new_node_id
end





function _refine!(volume::Volume{D},topo::Topology{D}) where D
    volume_id = volume.id

    if !is_root(volume)
        mid_node = find_single_intersec(get_volume_node_ids(topo),volume.childs)
        return mid_node
    end

    # @assert is_active(volume) "Inactive elements can't be refined. You are trying to refine volume $volume"
  
    volume.is_root = false
    volume_nodes = get_volume_node_ids(topo,volume_id)
    face_to_center = MutableSmallDict{20,Int,Int}()

    for area_id in get_volume_area_ids(topo,volume_id)
        area = get_areas(topo)[area_id]
        face_to_center[area_id] = _refine!(area,topo) 
    end
    
    was_once_refined = !isempty(volume.childs)
    if !was_once_refined
        # new_node_coords = mean(get_nodes(topo)[volume_nodes])
        new_node_coords = mean(get_nodes(topo)[n_id] for n_id in volume_nodes)
        new_node_id = add_node!(new_node_coords,topo)
    else
        new_node_id = find_single_intersec(get_volume_node_ids(topo),volume.childs)
        activate!(get_nodes(topo)[new_node_id])
        # reactivate childs and inner faces/edges
        for child_id in volume.childs
            activate!(get_volumes(topo)[child_id])
            for face_id in get_volume_area_ids(topo, child_id)
                activate!(get_areas(topo)[face_id])
                for edge_id in get_area_edge_ids(topo, face_id)
                    activate!(get_edges(topo)[edge_id])
                end
            end
        end
        return new_node_id
    end

    # Precompute face centers and edges
    face_ids = get_volume_area_ids(topo, volume_id)
    fid_to_corner_edge_ids = MutableSmallDict{20,Int,Tuple{Int,Int}}()
    boundary_quad_ids = @MVector zeros(Int,40)
    # For each corner node of the parent volume, create one child polyhedron
    for corner_vid in volume_nodes
        boundary_quad_counter = 1
        # Incident parent faces to this corner
        # For each incident face, find the two parent edges touching the corner and their midpoints
        for fid ∈ face_ids
            corner_vid ∈ get_area_node_ids(topo, fid) || continue

            quad_nodes = get_area_node_ids(topo, fid) 
            if corner_vid ∈ quad_nodes
                # push!(boundary_quads, quad_nodes)
                for child_face_id in get_areas(topo)[fid].childs 
                    child_quad_nodes = get_area_node_ids(topo, child_face_id)
                    if corner_vid ∈ child_quad_nodes
                        # push!(boundary_quads, child_quad_nodes)
                        boundary_quad_ids[boundary_quad_counter] = child_face_id
                        boundary_quad_counter += 1
                    end 
                end
            end 
        end

        # Internal faces: one quad per corner-edge using the two adjacent face centers
        # Collect all unique parent edges incident to the corner across incident faces

        empty!(fid_to_corner_edge_ids)

        for fid ∈ face_ids
            corner_vid ∈ get_area_node_ids(topo, fid) || continue
            for eid in get_area_edge_ids(topo, fid)
                if corner_vid in get_edge_node_ids(topo, eid)
                    edge_faces = get(fid_to_corner_edge_ids,eid,(-1,-1))
                    if edge_faces[1] == -1
                        fid_to_corner_edge_ids[eid] = (fid,edge_faces[2])
                    else
                        fid_to_corner_edge_ids[eid] = (edge_faces[1],fid)
                    end
                end
            end
        end
        for eid in keys(fid_to_corner_edge_ids)
            mid = find_single_intersec(get_edge_node_ids(topo),get_edges(topo)[eid].childs)
            
            adj_faces = fid_to_corner_edge_ids[eid]
            @assert length(adj_faces) == 2 "Corner edge must belong to exactly two incident faces."
            f1c = face_to_center[adj_faces[1]]
            f2c = face_to_center[adj_faces[2]]
            boundary_quad_ids[boundary_quad_counter] = add_area!(SA[mid, f1c, new_node_id, f2c],topo,0,volume.refinement_level)
            boundary_quad_counter += 1
            
        end

        nvolume_id = add_volume!(@view(boundary_quad_ids[1:boundary_quad_counter-1]), topo, volume_id, volume.refinement_level)
        push!(volume.childs, nvolume_id)
    end

    return new_node_id
end


