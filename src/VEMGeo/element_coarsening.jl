function downstream_deactivate!(mf::NManifold{N},topo::Topology{D}) where {N,D}

    deactivate!(mf)
    mf.is_root = true
    for child_id in mf.childs
        manifolds = get_geo_vec(topo,Val(N))
        !is_active(manifolds[child_id]) && continue 
        
        downstream_deactivate!(manifolds[child_id],topo) 
        mid_node_id = find_single_intersec(topo.connectivity[1,N],mf.childs)
        deactivate!(get_nodes(topo)[mid_node_id])
    end
end 


function _coarsen!(area::Area{D}, topo::Topology{D}) where D
    parent_id = area.parent_id
    @assert parent_id != 0 "Area has no parent"
    parent = get_areas(topo)[parent_id]
    !is_active(parent) && return

    bottom_up_coarsen_children!(parent, topo)
    mid_node_id = reactivate_parent_and_get_center!(parent, topo)
    deactivate_submanifolds_touching_center!(parent, topo, mid_node_id)
    for child_id in parent.childs
        deactivate!(get_areas(topo)[child_id])
    end
end



function _coarsen!(volume::Volume{D}, topo::Topology{D}) where D
    parent_id = volume.parent_id
    @assert parent_id != 0 "Volume has no parent"
    parent = get_volumes(topo)[parent_id]
    !is_active(parent) && return

    bottom_up_coarsen_children!(parent, topo)
    mid_node_id = reactivate_parent_and_get_center!(parent, topo)
    deactivate_submanifolds_touching_center!(parent, topo, mid_node_id)
    for child_id in parent.childs
        deactivate!(get_volumes(topo)[child_id])
    end
end


# ===================== Shared helpers =====================

function bottom_up_coarsen_children!(parent::NManifold{N}, topo::Topology{D}) where {N,D}
    manifolds = get_geo_vec(topo, Val(N))
    for child_id in parent.childs
        child = manifolds[child_id]
        is_root(child) && continue
        for sub_child_id in child.childs
            _coarsen!(manifolds[sub_child_id], topo)
        end
    end
    return nothing
end

function reactivate_parent_and_get_center!(parent::NManifold{N}, topo::Topology{D}) where {N,D}
    parent.is_root = true
    activate!(parent)
    mid_node_id = find_single_intersec(topo.connectivity[1, N], parent.childs)
    deactivate!(get_nodes(topo)[mid_node_id])
    return mid_node_id
end

function deactivate_submanifolds_touching_center!(parent::Area{D}, topo::Topology{D}, mid_node_id::Int) where D
    for child_id in parent.childs
        for edge_id in get_area_edge_ids(topo, child_id)
            n1, n2 = get_edge_node_ids(topo, edge_id)
            if n1 == mid_node_id || n2 == mid_node_id
                downstream_deactivate!(get_edges(topo)[edge_id], topo)
            end
        end
    end
    return nothing
end

function deactivate_submanifolds_touching_center!(
    parent::Volume{D}, 
    topo::Topology{D}, 
    mid_node_id::Int) where D

    for child_id in parent.childs
        for face_id in get_volume_area_ids(topo, child_id)
            if mid_node_id in get_area_node_ids(topo, face_id)
                downstream_deactivate!(get_areas(topo)[face_id], topo)
                for edge_id in get_area_edge_ids(topo, face_id)
                    n1id, n2id = get_edge_node_ids(topo, edge_id) 
                    if n1id == mid_node_id || n2id == mid_node_id
                        downstream_deactivate!(get_edges(topo)[edge_id], topo)
                    end
                end
            end
        end
    end
    return nothing
end