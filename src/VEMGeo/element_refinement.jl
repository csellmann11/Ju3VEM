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
    nodes = get_nodes(topo)
    edges = get_edges(topo)

    if !was_once_refined #standard case
        mid_node_coords = mean(nodes[n_id] for n_id in area_nodes)
        new_node_id     = add_node!(mid_node_coords,topo)

    else
        new_node_id = find_single_intersec(get_area_node_ids(topo),area.childs)
        activate!(nodes[new_node_id])
        # reactivate childs and inner edges
        for edge_id in get_area_edge_ids(topo,area_id)
            edge = edges[edge_id]
            _refine!(edge,topo)
        end

        for child_id in area.childs
            activate!(get_areas(topo)[child_id])
            for edge_id in get_area_edge_ids(topo,child_id)
                activate!(edges[edge_id])
            end
        end
        return new_node_id
    end

    @no_escape begin

        new_geo_nodes = @alloc(Int,2*length(area_nodes))
        
        for (edge_count,edgeid) in enumerate(get_area_edge_ids(topo,area_id))
            edge = edges[edgeid] 
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


function find_common_in_4(i1::Integer, i2::Integer, i3::Integer, i4::Integer)
    i1 == i2 && return i1
    i1 == i3 && return i1
    i1 == i4 && return i1
    i2 == i3 && return i2
    i2 == i4 && return i2
    i3 == i4 && return i3
    return -1
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
    full_area_node_ids = get_area_node_ids(topo)
    full_edge_node_ids = get_edge_node_ids(topo)
    areas              = get_areas(topo)
    for corner_vid in volume_nodes
        boundary_quad_counter = 1
        
        # Single pass: collect boundary child faces AND build edge-to-faces mapping
        empty!(fid_to_corner_edge_ids)
        
        for fid ∈ face_ids

            face_nodes = full_area_node_ids[fid]
            corner_vid ∈ face_nodes || continue
            
            # Collect boundary child faces containing this corner
            for child_face_id in areas[fid].childs 
                if corner_vid ∈ full_area_node_ids[child_face_id]
                    boundary_quad_ids[boundary_quad_counter] = child_face_id
                    boundary_quad_counter += 1
                end 
            end
            
            # Build edge-to-faces mapping for internal faces
            for eid in get_area_edge_ids(topo, fid)
                if corner_vid ∈ full_edge_node_ids[eid]
                    edge_faces = get(fid_to_corner_edge_ids, eid, (-1,-1))
                    fid_to_corner_edge_ids[eid] = edge_faces[1] == -1 ? (fid, edge_faces[2]) : (edge_faces[1], fid)
                end
            end
        end
        
        # Create internal faces connecting corner to volume center
        for eid in keys(fid_to_corner_edge_ids)
            edge = get_edges(topo)[eid]
            # mid = find_single_intersec(get_edge_node_ids(topo), get_edges(topo)[eid].childs)
            nc11,nc12 = full_edge_node_ids[edge.childs[1]]
            nc21,nc22 = full_edge_node_ids[edge.childs[2]]
            mid = find_common_in_4(nc11,nc12,nc21,nc22)
            # mid = find_single_intersec(nodes_child1, nodes_child2)
            adj_faces = fid_to_corner_edge_ids[eid]
            @assert length(adj_faces) == 2 "Corner edge must belong to exactly two incident faces."
            
            f1c = face_to_center[adj_faces[1]]
            f2c = face_to_center[adj_faces[2]]
            boundary_quad_ids[boundary_quad_counter] = add_area!(SA[mid, f1c, new_node_id, f2c], topo, 0, volume.refinement_level)
            
            boundary_quad_counter += 1
        end

        nvolume_id = add_volume!(@view(boundary_quad_ids[1:boundary_quad_counter-1]), topo, volume_id, volume.refinement_level)
        push!(volume.childs, nvolume_id)
    end

    return new_node_id
end





############################################################################################################
#################### Holistic Refinement ##################################################################
############################################################################################################

# function _refine!(area::Area{D},topo::Topology{D},edge_mid_node_ids::Vector{Int32}) where D
#     area_id = area.id

#     if !is_root(area)
#         mid_node_id = find_single_intersec(get_area_node_ids(topo),area.childs)
#         return mid_node_id
#     end

#     area.is_root = false
#     area_nodes = get_area_node_ids(topo,area_id)



 
    
#     was_once_refined = !isempty(area.childs)
    
#     if !was_once_refined #standard case
#         mid_node_coords = mean(get_nodes(topo)[n_id] for n_id in area_nodes)
#         new_node_id     = add_node!(mid_node_coords,topo)

#     else
#         new_node_id = find_single_intersec(get_area_node_ids(topo),area.childs)
#         activate!(get_nodes(topo)[new_node_id])
#         # reactivate childs and inner edges

#         for child_id in area.childs
#             activate!(get_areas(topo)[child_id])
#             for edge_id in get_area_edge_ids(topo,child_id)
#                 activate!(get_edges(topo)[edge_id])
#             end
#         end
#         return new_node_id
#     end

#     @no_escape begin

#         new_geo_nodes = @alloc(Int,2*length(area_nodes))
#         for (edge_count,edgeid) in enumerate(get_area_edge_ids(topo,area_id))
#             edge_mid_node_id = edge_mid_node_ids[edgeid]
#             idx = 2edge_count
#             new_geo_nodes[idx]   = edge_mid_node_id
#             new_geo_nodes[idx-1] = area_nodes[edge_count]
#         end

#         # Create quadrilaterals by connecting barycenter to refined edge nodes
#         sizehint!(area.childs, length(area_nodes))
#         first_node = new_geo_nodes[end]
#         for i in eachindex(area_nodes)
#             # For first sub-area: use last node, otherwise use node at 2i-2
#             # Create quadrilateral: [first_node, edge_node1, edge_node2, barycenter]
#             quadrilateral_nodes = SVector(first_node, new_geo_nodes[2i-1], new_geo_nodes[2i], new_node_id)

#             narea_id = add_area!(quadrilateral_nodes, topo, area_id, area.refinement_level)
#             push!(area.childs, narea_id)
#             first_node = new_geo_nodes[2i]
#         end
#         nothing
#     end
#     return new_node_id
# end


# function get_create_corner_faces!(
#     boundary_quad_ids::AbstractVector{<:Integer},
#     node_to_faces::Dict{Int32,Vector{Int32}},
#     node_to_edges::Dict{Int32,Vector{Int32}},
#     edge_mid_node_ids::Vector{Int32},
#     face_mid_node_ids::Vector{Int32},
#     topo::Topology{D},
#     volume::Volume{D},
#     el_center_node_id::Integer,
#     corner_vid::Integer,
#     volume_face_ids::AbstractVector{<:Integer},
#     volume_node_ids::AbstractVector{<:Integer}) where D
    
#     boundary_quad_counter = 1

#     corner_faces = node_to_faces[corner_vid]
#     for child_face_id in corner_faces
#         parent_id = get_areas(topo)[child_face_id].parent_id
#         parent_id in volume_face_ids || continue
#         boundary_quad_ids[boundary_quad_counter] = child_face_id
#         boundary_quad_counter += 1
#     end


#     adj_edge_ids = node_to_edges[corner_vid]
#     edge_nodes   = get_edge_node_ids(topo) 
#     areas        = get_areas(topo)
#     area_node_ids = get_area_node_ids(topo) 
#     for edge_id in adj_edge_ids 
#         _n1_id,_n2_id = edge_nodes[edge_id]

#         n2_id = _n1_id == corner_vid ? _n2_id : _n1_id
#         n2_id in volume_node_ids || continue
#         mid_node_id = edge_mid_node_ids[edge_id]
#         f1c = -1; f2c = -1;
#         for child_face_id in corner_faces
#             child_node_ids = area_node_ids[child_face_id] 
#             mid_node_id ∈ child_node_ids || continue
#             parent_id = areas[child_face_id].parent_id 
#             fc = face_mid_node_ids[parent_id]
#             f1c == -1 ? f1c = fc : f2c = fc
#         end
#         boundary_quad_ids[boundary_quad_counter] = add_area!(SA[mid_node_id, f1c, el_center_node_id, f2c], topo, 0, volume.refinement_level)
#         boundary_quad_counter += 1
#     end
#     return boundary_quad_counter
# end


# function _refine!(volume::Volume{D},topo::Topology{D},
#     node_to_faces::Dict{Int32,Vector{Int32}},
#     node_to_edges::Dict{Int32,Vector{Int32}},
#     face_mid_node_ids::Vector{Int32},
#     edge_mid_node_ids::Vector{Int32}
#     ) where D
#     volume_id = volume.id

#     if !is_root(volume)
#         mid_node = find_single_intersec(get_volume_node_ids(topo),volume.childs)
#         return mid_node
#     end

#     # @assert is_active(volume) "Inactive elements can't be refined. You are trying to refine volume $volume"
  
#     volume.is_root = false
#     volume_nodes   = get_volume_node_ids(topo,volume_id)

#     new_node_id = -1
#     if isempty(volume.childs) # Volume was never refined
#         new_node_coords = mean(get_nodes(topo)[n_id] for n_id in volume_nodes)
#         new_node_id = add_node!(new_node_coords,topo)
#     else
#         new_node_id = find_single_intersec(get_volume_node_ids(topo),volume.childs)
#         activate!(get_nodes(topo)[new_node_id])

#         # reactivate childs and inner faces/edges #TODO: Check why this is different then for areas
#         for child_id in volume.childs
#             activate!(get_volumes(topo)[child_id])
#             for face_id in get_volume_area_ids(topo, child_id)
#                 activate!(get_areas(topo)[face_id])
#                 for edge_id in get_area_edge_ids(topo, face_id)
#                     activate!(get_edges(topo)[edge_id])
#                 end
#             end
#         end 
#         return new_node_id
#     end

#     # For each corner node of the parent volume, create one child polyhedron
#     # fid_to_corner_edge_ids = MutableSmallDict{20,Int,Tuple{Int,Int}}()
#     boundary_quad_ids = @MVector zeros(Int,40)
#     volume_face_ids = get_volume_area_ids(topo, volume_id)
#     sizehint!(volume.childs, length(volume_nodes))
  
#     for corner_vid in volume_nodes
#         boundary_quad_counter = get_create_corner_faces!(
#             boundary_quad_ids,
#             node_to_faces,
#             node_to_edges,
#             edge_mid_node_ids,
#             face_mid_node_ids,
#             topo,
#             volume,
#             new_node_id,
#             corner_vid,
#             volume_face_ids,
#             volume_nodes)

#         nvolume_id = add_volume!(boundary_quad_ids[1:boundary_quad_counter-1], topo, volume_id, volume.refinement_level)
#         push!(volume.childs, nvolume_id)
#     end

#     return new_node_id
# end




# function refine!(topo::Topology{3}, elements_to_refine::BitVector) 

#     node_to_faces   = Dict{Int32,Vector{Int32}}()  #
#     node_to_edges  = Dict{Int32,Vector{Int32}}()  # this needs to be node_to_edges

#     edges_to_refine = BitVector(false for _ in get_edges(topo))  #
#     faces_to_refine = BitVector(false for _ in get_areas(topo))  #

#     edge_mid_node_ids = Vector{Int32}(undef, length(get_edges(topo)))  #
#     face_mid_node_ids = Vector{Int32}(undef, length(get_areas(topo)))  #

#     for element in RootIterator{4}(topo)
#         el_id = element.id 
#         elements_to_refine[el_id] || continue  

  
#         face_ids = get_volume_area_ids(topo, el_id)

#         for face_id in face_ids
#             faces_to_refine[face_id] = true

#             node_ids = get_area_node_ids(topo, face_id)
#             for (i,n_id) in enumerate(node_ids)
#                 ip1 = get_next_idx(node_ids,i)
#                 np1_id = node_ids[ip1]
#                 edge_id = get_edge(n_id,np1_id,topo) |> get_id 
#                 edges_to_refine[edge_id] = true
#             end
#         end
#     end


#     @time "edge_refine_loop" for i in eachindex(edges_to_refine)
#         edges_to_refine[i] || continue
        
#         n1_id, n2_id = get_edge_node_ids(topo, i)
#         push!(get!(() -> Int32[], node_to_edges, n1_id), i)
#         push!(get!(() -> Int32[], node_to_edges, n2_id), i)
        
#         mid_node_id = _refine!(get_edges(topo)[i],topo)
#         edge_mid_node_ids[i] = mid_node_id
#     end

#     @time "face_refine_loop" for i in eachindex(faces_to_refine)
#         faces_to_refine[i] || continue

#         face = get_areas(topo)[i]
#         mid_node_id = _refine!(face,topo,edge_mid_node_ids)
#         face_mid_node_ids[i] = mid_node_id
#         for (j,node_id) in enumerate(get_area_node_ids(topo, i))
#             child_j = face.childs[j]
#             push!(get!(() -> Int32[], node_to_faces, node_id), child_j)
#         end
#     end


#     @time "el_refine_loop" for i in eachindex(elements_to_refine)
#         elements_to_refine[i] || continue

#         el = get_volumes(topo)[i]
#        _refine!(el,topo,node_to_faces,node_to_edges,face_mid_node_ids,edge_mid_node_ids)
      
#     #    break


#     end


#     # el = get_volumes(topo)[1]
#     #  _refine!(el,topo,node_to_faces,node_to_edges,area_to_element,face_mid_node_ids,edge_mid_node_ids)



# end