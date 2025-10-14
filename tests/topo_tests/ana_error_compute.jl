import Ferrite 
import Ferrite.Tensors as TS
function compute_error(u_ana::F,
    cv::CellValues{D,U,ET},
    u::AbstractVector{Float64}) where {D,U,F<:Function,K,ET<:ElType{K}}

    base3d = get_base(BaseInfo{3,K,U}())

    ip = Ferrite.Lagrange{Ferrite.RefTetrahedron,1}()^U 
    qr = Ferrite.QuadratureRule{Ferrite.RefTetrahedron}(5)
    fr_cv = Ferrite.CellValues(qr,ip)


    topo = cv.mesh.topo

    l2_error = 0.0
    # h1_error = 0.0

    first = true

    for element in RootIterator{4}(cv.mesh.topo)
        reinit!(element.id,cv)

        hvol = cv.volume_data.hvol
        bc = cv.volume_data.vol_bc
        abs_volume = cv.volume_data.integrals[1]
        vol_data = cv.volume_data

        
        proj_s, proj = create_volume_vem_projectors(
            element.id,cv.mesh,cv.volume_data,cv.facedata_col,cv.vnm)
 

        dofs = Ju3VEM.VEMUtils.get_cell_dofs(cv)
        uel  = u[dofs] 

        nm = cv.vnm 
        for (node_id,local_id) in nm.map
            uel[local_id] = u[node_id]
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
                quad_coord_scaled = (quad_coord .- bc)/hvol
                
                uπx   = uπ(quad_coord_scaled)
                l2_error += norm(u_anax .- uπx)^2 * dΩ
            end

        end


    end

    return sqrt(l2_error)
end



function check_mesh(cv::CellValues{3,U}; is_dirichlet_boundary::F) where {U,F<:Function}

    mesh = cv.mesh
    topo = mesh.topo


    min_hvol = Inf
    max_hvol = 0.0 

    min_edge_length = Inf
    max_edge_length = 0.0

    max_coplanar_frac = 0.0

    max_face_nodes   = 0 
    max_vol_nodes    = 0 
    max_vol_faces    = 0
    all_faces_planar = 0

    node_id_to_element = Dict{Int,Vector{Int}}()
    
    for element in RootIterator{4}(topo)
        reinit!(element.id,cv)
        hvol = cv.volume_data.hvol
        bc = cv.volume_data.vol_bc

        n_nodes = Ju3VEM.VEMUtils.get_n_nodes(cv.vnm)
        max_vol_nodes = max(max_vol_nodes,n_nodes)

        vol_face_counter = Ref(0)
        max_face_nodes_vol = Ref(0)
        all_faces_planar_local = Ref(0)
        num_coplanar_faces_local = Ref(0)

        min_hvol = min(min_hvol,hvol)
        max_hvol = max(max_hvol,hvol)

        face_normals = SVector{3,Float64}[]

        nodes_seen_for_element = Set{Int}()

        Ju3VEM.VEMGeo.iterate_volume_areas(cv.facedata_col,topo,element.id) do _, face_data, _ 

            face_nodes_ids = face_data.face_node_ids.v.args[1]    
            face_nodes     = getindex.(Scalar(mesh.nodes),face_nodes_ids)

            for node_id in face_nodes_ids
                node_id ∈ nodes_seen_for_element && continue
                push!(nodes_seen_for_element,node_id)
                push!(get!(node_id_to_element,node_id,Int[]),element.id)
            end
            
            #Info: Check planarity 
            plane = face_data.dΩ.plane
            face_nodes2d = Ju3VEM.VEMGeo.project_to_2d_abs.(face_nodes,Scalar(plane)) 
            face_nodes3d_back = Ju3VEM.VEMGeo.project_to_3d.(face_nodes2d,Scalar(plane))
            if norm(face_nodes3d_back - face_nodes) > 1e-6 
                all_faces_planar_local[] += 1
            end

            #Info: Check if faces are coplanar 
            face_normal = plane.n
            for normal in face_normals 
                if norm(normal - face_normal) < 1e-6 || 
                    norm(normal - (-face_normal)) > 1e-6

                    num_coplanar_faces_local[] += 1
                end
            end

            for (i,node_id) in enumerate(face_nodes_ids)
                ip1 = get_next_idx(face_nodes_ids,i)
                edge_length = norm(face_nodes[i] - face_nodes[ip1])
                min_edge_length = min(min_edge_length,edge_length)
                max_edge_length = max(max_edge_length,edge_length)
            end

            push!(face_normals,face_normal)

            max_face_nodes_vol[] = max(max_face_nodes_vol[],length(face_nodes_ids))
            vol_face_counter[] += 1
        end

        frac = num_coplanar_faces_local[] / vol_face_counter[]
        max_coplanar_frac = max(max_coplanar_frac,frac)


        max_vol_faces = max(max_vol_faces,vol_face_counter[])
        max_face_nodes = max(max_face_nodes,max_face_nodes_vol[])
        all_faces_planar = all_faces_planar + all_faces_planar_local[]
    end

    for node_id in keys(node_id_to_element)
        node = mesh.nodes[node_id] 
        is_dirichlet_boundary(node) && continue 

        if length(node_id_to_element[node_id]) != 4 
            println("node $node_id has $(length(node_id_to_element[node_id])) elements")
            elements_unique = allunique(node_id_to_element[node_id])
            println("elements_unique: $elements_unique")
            display(node)
        end
    end



    println("finished checking mesh")
    println("="^50)
    println("General mesh information")
    println("="^50)
    println("The number of volumes is $(length(RootIterator{4}(topo)))")
    println("The number of faces is $(length(RootIterator{3}(topo)))")
    println("The number of nodes is $(length(get_nodes(topo)))")
    println("="^50)
    println("min_hvol: $min_hvol, max_hvol: $max_hvol")
    println("the maximum number of faces for a volume is $max_vol_faces")
    println("the maximum number of nodes for a face is $max_face_nodes")
    println("the maximum number of nodes for a volume is $max_vol_nodes")
    println("the number of faces which are not planar is $all_faces_planar")
    println("the maximum fraction of coplanar faces is $max_coplanar_frac")
    println("the minimum edge length is $min_edge_length, the maximum edge length is $max_edge_length")

    return nothing


end