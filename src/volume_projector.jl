using Ju3VEM

"""
    get_outward_normal(bcvol,face_data::FaceData{3,K,L}) where {K,L}

Returns the outward normal of the face w.r.t the volume center. 
!!!WARNING!!! Only works reliable for convex volumes.
"""
function get_outward_normal(bcvol,face_data::FaceData{3,K,L}) where {K,L}
    #TODO: move this function closer to geometric utils of the package
    plane = face_data.dΩ.plane
    normal = face_data.dΩ.n  

    face_bc3d = project_to_3d(face_data.bc,plane)
    face_to_vol_center = bcvol - face_bc3d 

    if dot(face_to_vol_center,normal) < 0
        return normal
    else
        return -normal
    end
end

function get_outward_normal(bcvol,face_data::FaceData{2,K,L}) where {K,L}
    return face_data.dΩ.plane.n
end




function create_volume_bmat(volume_id::Int,
    mesh::Mesh{3,ET},
    bcvol,hvol,
    facedata_col::Dict{Int,FaceData{3,K,L}},
    volume_node_mapping::ElementNodeMapping) where {K,L,ET<:ElType{K}}

    #TODO: write the face data structure
    topo = mesh.topo

    base_2d = get_base(BaseInfo{2,K,1}()).base
    base_3d = get_base(BaseInfo{3,K,1}()).base

    n_nodes = get_n_nodes(volume_node_mapping)
    bmat = zeros(length(base_3d),n_nodes)

    iterate_volume_areas(facedata_col,topo,volume_id) do root_area, face_data, topo

        face_normal      = get_outward_normal(bcvol,face_data)
        ΠsL2             = face_data.ΠsL2
        face_node_ids    = face_data.face_node_ids

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