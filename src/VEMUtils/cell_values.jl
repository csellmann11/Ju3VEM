

mutable struct CellValues{D,U,ET,L1,L2,V<:AbstractVector}

    const mesh::Mesh{D,ET}
    const vnm ::NodeID2LocalID

    const facedata_col::Dict{Int,FaceData{D,L1,V}}
    volume_data ::VolumeIntegralData{L2}

    const dh::DofHandler{D,U}
end



function _create_facedata_col(mesh::Mesh{D,StandardEl{K}}) where {D,K}
    K_MAX = max(2*K-1, 2)
    topo = mesh.topo


    FI_type = Core.Compiler.return_type(precompute_face_monomials,
            Tuple{Int,Topology{D},Val{K_MAX}})

    FD_Type = Core.Compiler.return_type(h1_projectors!,
            Tuple{Int,Mesh{3,StandardEl{K}},FI_type})
    facedata_col = Dict{Int,FD_Type}()

    for face in RootIterator{3}(topo)
        facedata_col[face.id] = h1_projectors!(
            face.id, mesh, 
            precompute_face_monomials(face.id, topo, Val(K_MAX))
            )
    end
    return facedata_col
end


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





struct ElementBaseFunctions{D,O,U,L,M<:StretchedMatrix}
    base::PolynomialBase{D,O,U,L}
    Π_star::M
end

function get_base_fun(ebf::ElementBaseFunctions, 
    idx::Int)
    return unit_sol_proj(ebf.base,idx,ebf.Π_star)
end

# iteration interface
function Base.iterate(ebf::ElementBaseFunctions, state::Int=1)
    N = size(ebf.Π_star, 2)
    state > N && return nothing
    return (get_base_fun(ebf, state), state + 1)
end

Base.length(ebf::ElementBaseFunctions) = size(ebf.Π_star, 2)

Base.eltype(::Type{ElementBaseFunctions{D,O,1,L,M}}) where {D,O,L,M} =
    Polynomial{Float64,D,L}  

Base.eltype(::Type{ElementBaseFunctions{D,O,U,L,M}}) where {D,O,U,L,M} =
    SVector{U,Polynomial{Float64,D,L}}