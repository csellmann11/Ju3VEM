using Ju3VEM 

using Ju3VEM.VEMGeo:FaceIntegralData

struct FaceData{D,K,L}
    face_node_ids::FlattenVecs{3,Int}
    dΩ           ::FaceIntegralData{D,K,L}
    ΠsL2         ::Matrix{Float64}
end
@inline Ju3VEM.VEMGeo.get_area(fd::FaceData) = get_area(fd.dΩ)
@inline Ju3VEM.VEMGeo.get_bc(fd::FaceData)   = get_bc(fd.dΩ)
@inline get_hf(fd::FaceData)   = fd.dΩ.hf

function create_B_mat(mesh::Mesh{D,ET},
    fd::FaceData) where {D,O,ET<:ElType{O}}

    area          = get_area(fd)
    hf            = get_hf(fd)
    center        = get_bc(fd)
    full_node_ids = fd.face_node_ids


    base     = get_base(BaseInfo{2,O,1}())
    mom_base = get_base(BaseInfo{2,O-2,1}())

    # full_node_ids   = get_iterative_area_node_ids(area,mesh)
    vertex_ids      = full_node_ids.v.args[1]
    edge_vertex_ids = full_node_ids.v.args[2]

    vertices      = @views mesh.nodes[vertex_ids]
    edge_vertices = @views mesh.nodes[edge_vertex_ids]


    num_vertices           = length(vertex_ids)
    num_edge_vertices      = length(edge_vertex_ids)


    num_nodes              = num_vertices + num_edge_vertices


    gw = get_gauss_lobatto(O+1)[2]

    B_mat = zeros(Float64,length(base),length(full_node_ids))
    @turbo B_mat .= 0.0
    if O == 1 
        B_mat[1,1:num_vertices] .= 1/num_vertices
    else
        B_mat[1,num_nodes+1] = 1
    end
    

    @inbounds for (idx,m) in enumerate(base)

        idx == 1 && continue
        
        # edge integral
        qp_count = 1
        # in the case of gradient values, m is already is gradient
        ∇m = ∇(m,hf)
        cnode_2d = SA[dot(vertices[1],fd.dΩ.u),dot(vertices[1],fd.dΩ.v)]
        for i in eachindex(vertices)
            ni = get_next_idx(vertices,i)

            # get current and next node
            nnode_2d = SA[dot(vertices[ni],fd.dΩ.u),dot(vertices[ni],fd.dΩ.v)]
            n, len = get_normal(nnode_2d,cnode_2d) 

            # influence of the mesh nodes. Each edge is toched by 2 nodes
            B_mat[idx,i]      += n ⋅ ∇m((cnode_2d-center)/hf)*len/2*gw[1]
            B_mat[idx,ni]     += n ⋅ ∇m((nnode_2d-center)/hf)*len/2*gw[end]

            # influence of the quadratur nodes.
            for j in 2:O
                quad_point = edge_vertices[qp_count]
                quad_point_2d = SA[dot(quad_point,fd.dΩ.u),dot(quad_point,fd.dΩ.v)]
                x = (quad_point_2d - center)/hf

                #INFO: For this integral the locations need to be lobatto points
                B_mat[idx,num_vertices+qp_count] = 
                        n ⋅ ∇m(x)*len/2*gw[j]
                qp_count += 1
            end
            cnode_2d = nnode_2d
        end

        # moment integral
        if sum(m.exp) < 2 continue end

        for d in 1:D
            ddm = ∂(m,hf,d,d)
            ddm.val == zero(ddm.val) && continue
            i = findfirst(m -> ddm.exp == m.exp,mom_base.base)
            B_mat[idx,num_nodes+i] -= ddm.val*area
        end
    end

    return B_mat
end


function create_D_mat(mesh::Mesh{D,ET},
    fd::FaceData) where {D,O,ET<:ElType{O}}

    base = get_base(BaseInfo{2,O,1}())

    area   = get_area(fd); 
    hf     = get_hf(fd); 
    center = get_bc(fd)

    vertex_ids      = fd.face_node_ids.v.args[1]
    edge_vertex_ids = fd.face_node_ids.v.args[2]
    moment_ids      = fd.face_node_ids.v.args[3]

    vertices        = @views mesh.nodes[vertex_ids]
    edge_vertices   = @views mesh.nodes[edge_vertex_ids]

    num_vertices = length(vertex_ids)
    num_edge_vertices = length(edge_vertex_ids)

    num_nodes = num_vertices + num_edge_vertices
    D_mat = zeros(Float64,num_nodes,length(base))
 

    for (idx,m) in enumerate(base)
        # Nodal evaluations
        for (i,node) in enumerate(Iterators.flatten((vertices,edge_vertices)))

            node_2d = SA[dot(node,fd.dΩ.u),dot(node,fd.dΩ.v)]
            D_mat[i,idx] = m(node_2d,center,hf)
        end

        for i in eachindex(moment_ids)
            D_mat[num_nodes+i,idx] = compute_face_integral_unshifted(m*base[i],fd.dΩ)/area
        end
    end
    return D_mat
end


#TODO: make face_id part of the face_data
function h1_projectors!(face_id::Int,mesh::Mesh{D,ET},
                       dΩ::FaceIntegralData) where {D,O,ET<:ElType{O}}

    topo = mesh.topo
    base   = get_base(BaseInfo{2,O,1}())

    full_node_ids = FlattenVecs{3,Int}()
    get_iterative_area_node_ids!(full_node_ids,get_areas(topo)[face_id],mesh)
    ΠsL2 = zeros(length(base),length(full_node_ids))
    face_data = FaceData(full_node_ids,dΩ,ΠsL2)


    D_mat  = create_D_mat(mesh,face_data)
    B_mat  = create_B_mat(mesh,face_data)
    

    G_mat  = SMatrix{length(base),length(base)}(B_mat * D_mat)
    invG_mat = inv(G_mat)


    Π_star = if O > 2
        zeros(size(B_mat))
    else 
        Π_star = face_data.ΠsL2 
    end

    Octavian.matmul!(Π_star,invG_mat,B_mat)

    if O > 2
        H_mat = create_H_mat(cv,b)
        C_mat = create_C_mat(cv,b,Π_star,H_mat) 
        Octavian.matmul!(Π_star,inv(H_mat),C_mat)
    end

    return face_data
end
