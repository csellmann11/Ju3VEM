import Ferrite 
import Ferrite.Tensors as TS
import Ferrite.ForwardDiff as FFD

function compute_error(u_ana::F,
    cv::CellValues{D,U,ET},
    u::AbstractVector{Float64}) where {D,U,F<:Function,K,ET<:ElType{K}}

    base3d = get_base(BaseInfo{3,K,U}())

    ip = Ferrite.Lagrange{Ferrite.RefTetrahedron,1}()^U 
    qr = Ferrite.QuadratureRule{Ferrite.RefTetrahedron}(5)
    fr_cv = Ferrite.CellValues(qr,ip)


    topo = cv.mesh.topo

    l2_error = 0.0
    h1_error = 0.0

    for element in RootIterator{4}(cv.mesh.topo)
        reinit!(element.id,cv)

        hvol = cv.volume_data.hvol
        bc = cv.volume_data.vol_bc 
        
        proj_s, _ = create_volume_vem_projectors(
            element.id,cv.mesh,cv.volume_data,cv.facedata_col,cv.vnm)
 

        dofs = get_cell_dofs(cv)
        uel  = u[dofs] 

        nm = cv.vnm 
        for (node_id,local_id) in nm.map
            dofs = get_dofs(cv.dh,node_id)
            for j in 1:U
                idx = (local_id-1)*U + j
                uel[idx] = u[dofs[j]]
            end
        end



        uπ   = sol_proj(base3d,uel,proj_s |> x -> stretch(x,Val(U)))

        tets_local, l2g = tetrahedralize_volume(topo, element.id)


        for i in eachindex(tets_local) 
            global_ids = l2g[tets_local[i]|> collect] 
            tet_nodes = TS.Vec{D}.(get_nodes(topo)[global_ids])
            Ferrite.reinit!(fr_cv,tet_nodes) 
            for qpoint in 1:Ferrite.getnquadpoints(fr_cv)
                dΩ = Ferrite.getdetJdV(fr_cv,qpoint) 
                quad_coord = Ferrite.spatial_coordinate(fr_cv,qpoint,tet_nodes)
                u_anax = u_ana(quad_coord)
                ∇u_anax = FFD.jacobian(x -> u_ana(x) |> SVector{U},quad_coord)
                quad_coord_scaled = (quad_coord .- bc)/hvol
                
                uπx   = uπ(quad_coord_scaled)
                ∇uπx = ∇x(uπ,hvol,quad_coord_scaled) |> SMatrix{U,D}
                l2_error += norm(u_anax .- uπx)^2 * dΩ
                h1_error += norm(∇u_anax .- ∇uπx)^2 * dΩ
            end

        end


    end

    return sqrt(l2_error), sqrt(h1_error)
end