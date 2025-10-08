mutable struct Assembler{T<:Real,D,U}
    
    const cv::CellValues{D,U}

    const max_entries::Int 
    const ndofs::Int
    count::Int

    const rows::Vector{Int}
    const cols::Vector{Int}
    const vals::Vector{T}
    const rhs::Vector{T}
end

#TODO: write this for 3d
function get_n_entries(cv::CellValues{2,U},order::Int) where U
    mesh = cv.mesh; topo = mesh.topo
    n_entries = 0 

    for element in RootIterator{2,3}(topo)
        nmn = length(topo.nids_col[element.id])
        nmoments = length(mesh.int_coords_connect[3][element.id])

        ndofs = (nmn + (order-1)*nmn + nmoments)*U
        n_entries += ndofs^2
    end
    n_entries
end

function get_n_entries(cv::CellValues{3,U},::Int) where U
    n_entries = 0
 
    nvol_nodes = Ref(0)

    for element in RootIterator{4}(cv.mesh.topo) 

        nvol_nodes[] = 0
        iterate_volume_areas(cv.facedata_col,cv.mesh.topo,element.id) do _, face_data, _
            n_face_nodes = length(face_data.face_node_ids)
            nvol_nodes[] += (n_face_nodes)
        end
        n_entries += (div(nvol_nodes[],2)*U)^2
        
    end
    n_entries[]
end

function Assembler{T}(cv::CellValues{D,U,ET},count = 1)  where {T<:Real,D,U,K,ET<:ElType{K}}
    max_entries = get_n_entries(cv,1)


    rows = Vector{Int}(undef, max_entries)
    cols = Vector{Int}(undef, max_entries)
    vals = Vector{T}(undef, max_entries)

    ndofs = length(cv.dh.dof_mapping)*U 
    rhs  = zeros(T, ndofs) 

    return Assembler(cv,max_entries, ndofs, count, rows, cols, vals, rhs)
end
Assembler(order, cv, count = 1) = Assembler{Float64}(order,cv,count)


@inline function local_ids_to_local_dofs(id::Integer,::DofHandler{D,U}) where {D,U}

    start_id = (id-1)*U + 1
    end_id   = id*U 
    return SVector{U}(start_id:end_id)
end

function local_assembly!(as::Assembler{T,D,U}, klocal::AbstractMatrix{T},
    flocal::AbstractVector{T}) where {T<:Real,D,U}
    
    @assert size(klocal,1) == size(flocal,1) == size(klocal,2) "Matrix and vector dimension mismatch"

    dh = as.cv.dh

    for (node_id_i,i) in as.cv.vnm.map
        dofs_i = get_dofs(dh, node_id_i) 
        local_dofs_i = local_ids_to_local_dofs(i, dh)
        as.rhs[dofs_i] += flocal[local_dofs_i]
         
        for (node_id_j,j) in as.cv.vnm.map
            dofs_j = get_dofs(dh, node_id_j)
            local_dofs_j = local_ids_to_local_dofs(j, dh)
            
            for u in 1:U, v in 1:U
                as.rows[as.count] = dofs_i[u]
                as.cols[as.count] = dofs_j[v]
                as.vals[as.count] = klocal[local_dofs_i[u], local_dofs_j[v]]
                as.count += 1
            end
        end
    end
end

function assemble!(as::Assembler)
    n_entries = as.count - 1
    return @views sparse(as.rows[1:n_entries], 
                        as.cols[1:n_entries], 
                        as.vals[1:n_entries], 
                        as.ndofs, as.ndofs), as.rhs
end