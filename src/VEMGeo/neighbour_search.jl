# """
#     create_element_neighbour_list(topo::Topology{D}) where D

# Function create for each active element a list of its neighbours. 
# If there is no neighbour on one side of the element, the list contains -face_id of this position
# List contains neighbours by edges and also neighbours which only share an edge.
# """
# _push_to_map!(map::Dict{Int,Vector{Int}}, key::Int, val::Int) = (push!(get!(map, key, Int[]), val); nothing)

# function create_element_neighbour_list(topo::Topology{3}; n_neighs_min::Integer=16)

#     # Collect active root volumes
#     # active_vol_ids = Int[]
#     # for vol in RootIterator{4}(topo)
#     #     push!(active_vol_ids, vol.id)
#     # end

#     # Map face -> active volumes, edge -> active volumes, node -> active volumes
#     face_to_vols = Dict{Int,Vector{Int}}()
#     edge_to_vols = Dict{Int,Vector{Int}}()
#     node_to_vols = Dict{Int,Vector{Int}}()

#     second_ring_seen = Set{Int}()
#     for vol in RootIterator{4}(topo)
#         vid = vol.id
#         for face_id in get_volume_area_ids(topo, vid)
#             _push_to_map!(face_to_vols, face_id, vid)
#         end
#         seen_edges = Set{Int}()
#         for face_id in get_volume_area_ids(topo, vid)
#             for e in get_area_edge_ids(topo, face_id)
#                 if e in seen_edges
#                     continue
#                 end
#                 _push_to_map!(edge_to_vols, e, vid)
#                 push!(seen_edges, e)
#             end
#         end
#         # Build node -> volumes mapping
#         seen_nodes = Set{Int}()
#         for face_id in get_volume_area_ids(topo, vid)
#             for n in get_area_node_ids(topo, face_id)
#                 if n in seen_nodes
#                     continue
#                 end
#                 _push_to_map!(node_to_vols, n, vid)
#                 push!(seen_nodes, n)
#             end
#         end
#     end

#     # Build neighbor vectors as FlattenVecs per volume: [faces, edges, second-ring, node-only]
#     neighs_dict = Dict{Int,FlattenVecs{3,Int}}()
#     for vol in RootIterator{4}(topo)
#         vid = vol.id
#         faces = get_volume_area_ids(topo, vid)

#         # Face neighbors 
#         face_neighs = Int[]
#         for face_id in faces
#             vols = get(face_to_vols, face_id, Int[])
#             if length(vols) >= 2
#                 for other in vols
#                     if other != vid
#                         push!(face_neighs, other)
#                         break
#                     end
#                 end
#             else
#                 push!(face_neighs, -face_id)
#             end
            
#         end
#         unique!(face_neighs)

#         # Edge-only neighbors (exclude already added face neighbors)
#         face_set = Set{Int}(face_neighs)
#         edge_neighs = Int[]
#         seen_edges = Set{Int}()
#         for face_id in faces
#             for e in get_area_edge_ids(topo, face_id)
#                 if e in seen_edges
#                     continue
#                 end
#                 push!(seen_edges, e)
#                 for other in get(edge_to_vols, e, Int[])
#                     if other != vid && !(other in face_set)
#                         push!(edge_neighs, other)
#                     end
#                 end
#             end
#         end
#         unique!(edge_neighs)
#         neighs_dict[vid] = FlattenVecs(face_neighs, edge_neighs,Int[])
#     end


#     for (vid,neighs) in neighs_dict
#         empty!(second_ring_seen)
#         union!(second_ring_seen, neighs)
#         second_ring = neighs.v.args[3]
#         if length(neighs) < n_neighs_min
#             for neight_i in neighs.v.args[1]

#                 if neight_i < 0 
#                     push!(second_ring, neight_i)
#                     continue 
#                 end

#                 neighs_neighs = neighs_dict[neight_i].v.args[1]
#                 for neigh_j in neighs_neighs
#                     neigh_j == vid && continue  
#                     neigh_j ∈ second_ring_seen && continue 
#                     push!(second_ring_seen, neigh_j)
#                     push!(second_ring, neigh_j)
#                 end
#             end
#         end

#     end


#     return neighs_dict,face_to_vols
# end

# function create_element_neighbour_list_fast(topo::Topology{3}; n_neighs_min::Int=16)
#     # Active volume ids and a compact index map (dense only for volumes)

#     active_vol_ids = [vol.id for vol in RootIterator{4}(topo)]

#     vid_to_idx = Dict{Int,Int}()
#     for (i, vid) in enumerate(active_vol_ids)
#         vid_to_idx[vid] = i
#     end

#     # Sparse maps: face -> active vols, edge -> active vols
#     face_to_vols = Dict{Int,Vector{Int}}()
#     edge_to_vols = Dict{Int,Vector{Int}}()

#     # Global edge-stamp dict to avoid per-volume Set allocations
#     edge_seen_stamp = Dict{Int,Int}()  # edge_id -> last_stamp_seen
#     edge_stamp = 0

#     # Pass A: build face_to_vols and edge_to_vols (dedup edges within a volume via stamps)
#     for vid in active_vol_ids
#         faces = get_volume_area_ids(topo, vid)
#         for f in faces
#             _push_to_map!(face_to_vols, f, vid)
#         end
#         edge_stamp += 1
#         s = edge_stamp
#         for f in faces
#             for e in get_area_edge_ids(topo, f)
#                 if get(edge_seen_stamp, e, 0) != s
#                     _push_to_map!(edge_to_vols, e, vid)
#                     edge_seen_stamp[e] = s
#                 end
#             end
#         end
#     end

#     # Stamps over compact volume indices: O(1) dedup without hashing
#     marks_included = zeros(Int, length(active_vol_ids))  # self + face included
#     marks_dedup    = zeros(Int, length(active_vol_ids))  # edge-only dedup
#     included_stamp = 0
#     dedup_stamp    = 0

#     # Pass B: face and edge neighbors (no Sets/unique!)
#     neighs_dict = Dict{Int,FlattenVecs{3,Int}}()
#     for vid in active_vol_ids
#         faces = get_volume_area_ids(topo, vid)

#         included_stamp += 1
#         marks_included[vid_to_idx[vid]] = included_stamp

#         face_neighs = Int[]
#         sizehint!(face_neighs, length(faces))

#         for f in faces
#             arr = get(face_to_vols, f, Int[])
#             if length(arr) >= 2
#                 nb = 0
#                 @inbounds for other in arr
#                     if other != vid
#                         nb = other
#                         break
#                     end
#                 end
#                 if nb != 0
#                     push!(face_neighs, nb)
#                     marks_included[vid_to_idx[nb]] = included_stamp
#                 else
#                     push!(face_neighs, -f)  # no "other" found
#                 end
#             else
#                 push!(face_neighs, -f)      # boundary face
#             end
#         end

#         edge_neighs = Int[]
#         dedup_stamp += 1
#         for f in faces
#             for e in get_area_edge_ids(topo, f)
#                 for nb in get(edge_to_vols, e, Int[])
#                     if nb != vid
#                         j = vid_to_idx[nb]
#                         if marks_included[j] != included_stamp && marks_dedup[j] != dedup_stamp
#                             marks_dedup[j] = dedup_stamp
#                             push!(edge_neighs, nb)
#                         end
#                     end
#                 end
#             end
#         end

#         neighs_dict[vid] = FlattenVecs(face_neighs, edge_neighs, Int[])
#     end

#     # Pass C: second ring from face-neighbors only (carry negative faces)
#     for (vid, packs) in neighs_dict
#         face_neighs = packs.v.args[1]
#         edge_neighs = packs.v.args[2]
#         ring2       = packs.v.args[3]

#         if length(face_neighs) + length(edge_neighs) < n_neighs_min
#             included_stamp += 1
#             marks_included[vid_to_idx[vid]] = included_stamp
#             for nb in face_neighs
#                 nb > 0 && (marks_included[vid_to_idx[nb]] = included_stamp)
#             end
#             for nb in edge_neighs
#                 marks_included[vid_to_idx[nb]] = included_stamp
#             end

#             for nb in face_neighs
#                 if nb < 0
#                     push!(ring2, nb)
#                     continue
#                 end
#                 nb_face = neighs_dict[nb].v.args[1]
#                 for nb2 in nb_face
#                     if nb2 > 0 && nb2 != vid
#                         j = vid_to_idx[nb2]
#                         if marks_included[j] != included_stamp
#                             marks_included[j] = included_stamp
#                             push!(ring2, nb2)
#                         end
#                     end
#                 end
#             end
#         end
#     end

#     return neighs_dict, face_to_vols
# end


# function create_element_neighbour_list(topo::Topology{2})
#     # Collect active root areas
#     active_area_ids = Int[]
#     for area in RootIterator{3}(topo)
#         push!(active_area_ids, area.id)
#     end

#     # Map edge -> active areas (for edge adjacency) and node -> active areas (vertex-only)
#     edge_to_areas = Dict{Int,Vector{Int}}()
#     node_to_areas = Dict{Int,Vector{Int}}()
#     for aid in active_area_ids
#         for e in get_area_edge_ids(topo, aid)
#             _push_to_map!(edge_to_areas, e, aid)
#         end
#         for n in get_area_node_ids(topo, aid)
#             _push_to_map!(node_to_areas, n, aid)
#         end
#     end

#     empty!(topo.el_neighs)
#     for aid in active_area_ids
#         edges = get_area_edge_ids(topo, aid)
#         neighs = Vector{Int}(undef, length(edges))

#         # Edge neighbors (or negative edge id if boundary)
#         for (i, e) in enumerate(edges)
#             arr = get(edge_to_areas, e, Int[])
#             nb = 0
#             if length(arr) >= 2
#                 for other in arr
#                     if other != aid
#                         nb = other
#                         break
#                     end
#                 end
#             end
#             neighs[i] = nb == 0 ? -e : nb
#         end

#         # Vertex-only neighbors (exclude already added edge neighbors)
#         edge_nb = Set{Int}(x for x in neighs if x > 0)
#         extras = Int[]
#         seen_nodes = Set{Int}()
#         for n in get_area_node_ids(topo, aid)
#             if n in seen_nodes
#                 continue
#             end
#             push!(seen_nodes, n)
#             for other in get(node_to_areas, n, Int[])
#                 if other != aid && !(other in edge_nb)
#                     push!(extras, other)
#                 end
#             end
#         end
#         if !isempty(extras)
#             unique!(extras)
#             neighs = vcat(neighs, extras)
#         end

#         topo.el_neighs[aid] = neighs
#     end
#     return topo.el_neighs
# end