

mutable struct CellValues{D,U,ET,L1,L2,V<:AbstractVector}

    const mesh::Mesh{D,ET}
    const vnm ::NodeID2LocalID

    const facedata_col::Dict{Int,FaceData{D,L1,V}}
    volume_data ::VolumeIntegralData{L2}

    const dh::DofHandler{D,U}
end



function _create_facedata_col(mesh::Mesh{D,StandardEl{K}}) where {D,K}
    K_MAX = max(2*K-1, 2)

    base = get_base(BaseInfo{2, K_MAX, 1}()).base
    topo = mesh.topo

    facedata_col = Dict(
        face.id => h1_projectors!(
            face.id, mesh, 
            precompute_face_monomials(face.id, topo, Val(K_MAX))
            ) for face in RootIterator{3}(topo)
    )

    #INFO: debugging util 
    for (id,fd) in facedata_col
        area = get_area(fd)
        # @assert area > 0 "area is zero for face $id with nodes $(fd.face_node_ids)"
        if area <= 0
            node_ids = fd.face_node_ids 
            nodes = get_nodes(topo)[node_ids]
            @show node_ids, nodes
            @show area
        end
    end

    # # Info: debugging check if face_plan 
    # for (id,fd) in facedata_col
    #     plane = fd.dΩ.plane 
    #     node_ids = fd.face_node_ids 
    #     nodes = get_nodes(topo)[node_ids]
        

    #     nodes_2d = project_to_2d_abs.(nodes,Ref(plane))
    #     nodes_recomp = project_to_3d.(nodes_2d,Ref(plane))
    #     for (node1,node2) in zip(nodes,nodes_recomp)
    #         if node1 ≉ node2
    #             display(node_ids)
    #             display(nodes)
    #             display(nodes_recomp)
    #             throw("face is not planar")
    #         end
    #     end
    # end


    return facedata_col
end


# function CellValues(
#     mesh::Mesh{D,ET},
#     vnm::NodeID2LocalID,
#     facedata_col::Dict{Int,FaceData{D,L1,V}},
#     volume_data::VolumeIntegralData{L2},
#     dh::DofHandler{D,U}) where {D,U,ET,L1,L2,V<:AbstractVector}

#     return CellValues{D,U,ET,L1,L2,V}(
#         mesh,vnm,facedata_col,volume_data,dh)
# end

function CellValues{U}(
    mesh::Mesh{D,ET}) where {D,U,K,ET<:ElType{K}}

    dh = DofHandler{U}(mesh)

    vnm = NodeID2LocalID(;sizehint = 30) 

    # base is just a dummy 3d base, even if problem is 2d
    base = get_base(BaseInfo{3,K,1}()).base
    volume_data = VolumeIntegralData{length(base)}(
        0.0,
        zero(SVector{3,Float64}),
        zero(SVector{length(base),Float64}))

    facedata_col = _create_facedata_col(mesh)
    V = typeof(first(values(facedata_col)).face_node_ids)

    return CellValues(mesh,vnm,facedata_col,volume_data,dh)
    # return CellValues{D,U,ET,length(base),length(base),V}(
    #     mesh,vnm,facedata_col,volume_data,dh)
end

function get_element_bc(cv::CellValues)
    return cv.volume_data.vol_bc
end

function get_element_diameter(cv::CellValues)
    return cv.volume_data.hvol
end

function get_element_area(cv::CellValues)
    return cv.volume_data.integrals[1]
end


function reinit!(el_id::Int,cv::CellValues{3})
    
    mesh = cv.mesh; topo = mesh.topo
    K = get_order(cv.mesh)

    vol_data = precompute_volume_monomials(
        el_id, topo, cv.facedata_col, Val(K))

    create_node_mapping!(cv.vnm,el_id,mesh,cv.facedata_col)

    cv.volume_data = vol_data

    return cv
end

@inline get_n_cell_dofs(cv::CellValues{D,U}) where {D,U} = U*get_n_nodes(cv.vnm)


