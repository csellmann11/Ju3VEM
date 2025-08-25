using Ju3VEM 

using Ju3VEM.VEMGeo:FaceIntegralData
using Ju3VEM.VEMGeo: D2FaceParametrization, project_to_2d_abs, project_to_2d_rel, project_to_3d

struct FaceData{D,K,L,FV<:FlattenVecs{3,Int}}
    face_node_ids::FV
    dΩ           ::FaceIntegralData{D,K,L}
    ΠsL2         ::FixedSizeMatrixDefault{Float64}
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

    plane = fd.dΩ.plane

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

    # B_mat = zeros(Float64,length(base),length(full_node_ids))
    B_mat = FixedSizeMatrix{Float64}(undef,length(base),length(full_node_ids))
    B_mat .= 0.0

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
        cnode_2d = project_to_2d_abs(vertices[1],plane)
        for i in eachindex(vertices)
            ni = get_next_idx(vertices,i)

            # get current and next node
            nnode_2d = project_to_2d_abs(vertices[ni],plane)
            n, len = get_normal(nnode_2d,cnode_2d) 

            # influence of the mesh nodes. Each edge is toched by 2 nodes
            B_mat[idx,i]      += n ⋅ ∇m((cnode_2d-center)/hf)*len/2*gw[1]
            B_mat[idx,ni]     += n ⋅ ∇m((nnode_2d-center)/hf)*len/2*gw[end]

            # influence of the quadratur nodes.
            for j in 2:O
                quad_point = edge_vertices[qp_count]
                quad_point_2d = project_to_2d_abs(quad_point,plane)
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
    plane = fd.dΩ.plane

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
    # D_mat = zeros(Float64,num_nodes,length(base))
    D_mat = FixedSizeMatrix{Float64}(undef,num_nodes,length(base))
 

    for (idx,m) in enumerate(base)
        # Nodal evaluations
        for (i,node) in enumerate(Iterators.flatten((vertices,edge_vertices)))

            node_2d = project_to_2d_abs(node,plane)
            D_mat[i,idx] = m(node_2d,center,hf)
        end

        for i in eachindex(moment_ids)
            D_mat[num_nodes+i,idx] = compute_face_integral_unshifted(m*base[i],fd.dΩ)/area
        end
    end
    return D_mat
end


#TODO: finish this function
function create_G_mat(mesh::Mesh{D,ET},
    fd::FaceData) where {D,O,ET<:ElType{O}}

    throw("not implemented")

    dΩ = fd.dΩ
    base = get_base(BaseInfo{2,O,1}())
    area = get_area(fd)
    center = get_bc(fd)
    hf = get_hf(fd)


    G_mat = MMatrix{length(base),length(base)}(undef)

    for (i,mi) in enumerate(base)
        i == 1 && continue
        for (j,mj) in enumerate(base)
            j == 1 && continue
            G_mat[i,j] = compute_face_integral_unshifted(mi*mj,dΩ)/area
        end
    end

    vertex_ids = fd.face_node_ids.v.args[1]
    if O == 1 
        for i in eachindex(vertex_ids)
            vertex = mesh.nodes[vertex_ids[i]]
            vertex_2d = project_to_2d_abs(vertex,dΩ.plane)
            for (j,mj) in enumerate(base)
                G_mat[1,j] += mj(vertex_2d,center,hf)/length(vertex_ids)
            end
        end
    else 


    end



    return G_mat

end
# @inline function static_matmul(A::AbstractMatrix{T}, 
#     B::AbstractMatrix{T}, ::Val{MN}) where {T,MN}
#     M,N = MN
#     @inbounds begin
#         return SMatrix{M,N,T}(
#             begin 
#                 s = zero(T)
#                 for k in axes(A,2)
#                     s += A[i,k] * B[k,j]
#                 end
#                 s
#             end
#             for i in 1:M, j in 1:N)
#     end
# end

@inline function static_matmul(A::AbstractMatrix{T}, B::AbstractMatrix{T}, ::Val{MN}) where {T,MN}
    M, N = MN
    C = MMatrix{M,N,T}(undef)
    @turbo for i in 1:M, j in 1:N
        s = zero(T)
        for k in axes(A,2)
            s += A[i,k] * B[k,j]
        end
        C[i,j] = s
    end
    return SMatrix(C)
end

#TODO: make face_id part of the face_data
function h1_projectors!(face_id::Int,mesh::Mesh{D,ET},
                       dΩ::FaceIntegralData) where {D,O,ET<:ElType{O}}

    topo = mesh.topo
    base   = get_base(BaseInfo{2,O,1}())

    face  = get_areas(topo)[face_id]
    full_node_ids = FlattenVecs{3,Int}()
    get_iterative_area_node_ids!(full_node_ids,face,mesh)

    # ΠsL2 = zeros(length(base),length(full_node_ids))
    ΠsL2 = FixedSizeMatrix{Float64}(undef,length(base),length(full_node_ids))
    ΠsL2 .= 0.0

    face_data = FaceData(full_node_ids,dΩ,ΠsL2)


    D_mat  = create_D_mat(mesh,face_data)
    B_mat  = create_B_mat(mesh,face_data)
    

    invG_mat = static_matmul(B_mat,D_mat,Val((length(base),length(base)))) |> inv


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
