using Ju3VEM

function create_volume_bmat(mesh::Mesh{3,ET},
    bcvol,hvol,
    facedata_col,
    volume_node_mapping::ElementNodeMapping) where {K,ET<:ElType{K}}

    #TODO: write the face data structure
    topo = mesh.topo

    base_2d = get_base(BaseInfo{2,K,1}()).base
    base_3d = get_base(BaseInfo{3,K,1}()).base

    n_nodes = get_n_nodes(volume_node_mapping)
    bmat = zeros(length(base_3d),n_nodes)

    iterate_volume_areas(topo,volume_id) do root_area, parent_area, volume_id, topo


        face_data        = facedata_col[root_area.id]
        face_normal      = get_outward_normal(bcvol,face_data)

        #TODO: write the face projectors
        ΠsL2             = face_data.ΠsL2

        #! maby collect those already when computing the face projector
        face_node_ids    = get_iterative_area_node_ids(root_area,mesh)

        integrated_vals = MVector{length(base_2d),Float64}(undef)
        for (i,m3d) in enumerate(base_3d)
            ∇m3d = ∇(m3d,hvol)
            integrated_vals .= 0.0
            for (j,m2d) in enumerate(base_2d)
                for d in 1:3
                    ∇m3d_nd = ∇m3d[d]*face_normal[d]
                    #TODO: write this function
                    integrated_vals[j] += get_face_integral(m2d,∇m3d_nd,face_data)
                end
            end
            for nc in eachindex(face_node_ids)
                val = @views dot(ΠsL2[nc,:],integrated_vals)
                node_id = face_node_ids[nc] 
                col_id  = volume_node_mapping[node_id]
                bmat[i,col_id] += val
            end
        end
    end
    return bmat
end