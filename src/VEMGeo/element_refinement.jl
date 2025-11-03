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



 
    
    was_once_refined = !isempty(area.childs)
    
    if !was_once_refined #standard case
        mid_node_coords = mean(get_nodes(topo)[n_id] for n_id in area_nodes)
        new_node_id     = add_node!(mid_node_coords,topo)

    else
        new_node_id = find_single_intersec(get_area_node_ids(topo),area.childs)
        activate!(get_nodes(topo)[new_node_id])
        # reactivate childs and inner edges
        for edge_id in get_area_edge_ids(topo,area_id)
            edge = get_edges(topo)[edge_id]
            _refine!(edge,topo)
        end

        for child_id in area.childs
            activate!(get_areas(topo)[child_id])
            for edge_id in get_area_edge_ids(topo,child_id)
                activate!(get_edges(topo)[edge_id])
            end
        end
        return new_node_id
    end

    @no_escape begin

        new_geo_nodes = @alloc(Int,2*length(area_nodes))
        for (edge_count,edgeid) in enumerate(get_area_edge_ids(topo,area_id))
            edge = get_edges(topo)[edgeid] 
            edge_mid_node_id = _refine!(edge,topo)
            idx = 2edge_count
            new_geo_nodes[idx]   = edge_mid_node_id
            new_geo_nodes[idx-1] = area_nodes[edge_count]
        end

        # Create quadrilaterals by connecting barycenter to refined edge nodes
        first_node = new_geo_nodes[end]
        for i in eachindex(area_nodes)
            # For first sub-area: use last node, otherwise use node at 2i-2
            # Create quadrilateral: [first_node, edge_node1, edge_node2, barycenter]
            quadrilateral_nodes = SVector(first_node, new_geo_nodes[2i-1], new_geo_nodes[2i], new_node_id)

            narea_id = add_area!(quadrilateral_nodes, topo, area_id, area.refinement_level)
            push!(area.childs, narea_id)
            first_node = new_geo_nodes[2i]
        end
        nothing
    end
    return new_node_id
end





function _refine!(volume::Volume{D},topo::Topology{D};
    non_planar = true) where D
    volume_id = volume.id

    if !is_root(volume)
        mid_node = find_single_intersec(get_volume_node_ids(topo),volume.childs)
        return mid_node
    end

    # @assert is_active(volume) "Inactive elements can't be refined. You are trying to refine volume $volume"
  
    volume.is_root = false
    volume_nodes = get_volume_node_ids(topo,volume_id)
    face_ids = get_volume_area_ids(topo, volume_id)
    face_to_center = MutableSmallDict{20,Int,Int}()

    # Refine all faces and store their centers
    for area_id in face_ids
        area = get_areas(topo)[area_id]
        face_to_center[area_id] = _refine!(area,topo) 
    end
 

    if isempty(volume.childs) # Volume was never refined
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

    # For each corner node of the parent volume, create one child polyhedron
    fid_to_corner_edge_ids = MutableSmallDict{20,Int,Tuple{Int,Int}}()
    boundary_quad_ids = @MVector zeros(Int,40)
    
    for corner_vid in volume_nodes
        boundary_quad_counter = 1
        
        # Single pass: collect boundary child faces AND build edge-to-faces mapping
        empty!(fid_to_corner_edge_ids)
        
        for fid ∈ face_ids
            face_nodes = get_area_node_ids(topo, fid)
            corner_vid ∈ face_nodes || continue
            
            # Collect boundary child faces containing this corner
            for child_face_id in get_areas(topo)[fid].childs 
                if corner_vid ∈ get_area_node_ids(topo, child_face_id)
                    boundary_quad_ids[boundary_quad_counter] = child_face_id
                    boundary_quad_counter += 1
                end 
            end
            
            # Build edge-to-faces mapping for internal faces
            for eid in get_area_edge_ids(topo, fid)
                if corner_vid ∈ get_edge_node_ids(topo, eid)
                    edge_faces = get(fid_to_corner_edge_ids, eid, (-1,-1))
                    fid_to_corner_edge_ids[eid] = edge_faces[1] == -1 ? (fid, edge_faces[2]) : (edge_faces[1], fid)
                end
            end
        end
        
        # Create internal faces connecting corner to volume center
        for eid in keys(fid_to_corner_edge_ids)
            edge = get_edges(topo)[eid]
            # mid = find_single_intersec(get_edge_node_ids(topo), get_edges(topo)[eid].childs)
            nodes_child1 = get_edge_node_ids(topo,edge.childs[1])
            nodes_child2 = get_edge_node_ids(topo,edge.childs[2])
            mid = find_single_intersec(nodes_child1, nodes_child2)
            adj_faces = fid_to_corner_edge_ids[eid]
            @assert length(adj_faces) == 2 "Corner edge must belong to exactly two incident faces."
            
            f1c = face_to_center[adj_faces[1]]
            f2c = face_to_center[adj_faces[2]]
            boundary_quad_ids[boundary_quad_counter] = add_area!(SA[mid, f1c, new_node_id, f2c], topo, 0, volume.refinement_level)
            
            # if non_planar
            #     last_area_id = boundary_quad_ids[boundary_quad_counter]
            #     last_area_added = get_areas(topo)[last_area_id]
 
            #     # triangulate the face if its non-planar
            #     tria1_id = add_area!(SA[mid, f1c, new_node_id], topo, last_area_id, volume.refinement_level)
            #     tria2_id = add_area!(SA[new_node_id, f2c, mid], topo, last_area_id, volume.refinement_level)
            #     append!(last_area_added.childs, [tria1_id, tria2_id])
            #     last_area_added.is_root = false

            #     boundary_quad_ids[boundary_quad_counter] = tria1_id
            #     boundary_quad_counter += 1
            #     boundary_quad_ids[boundary_quad_counter] = tria2_id
            # end
            boundary_quad_counter += 1
        end

        nvolume_id = add_volume!(@view(boundary_quad_ids[1:boundary_quad_counter-1]), topo, volume_id, volume.refinement_level)
        push!(volume.childs, nvolume_id)
    end

    return new_node_id
end

