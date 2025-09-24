# using WriteVTK
# using Ju3VEM.VEMGeo: get_iterative_area_vertex_ids,get_unique_values
# # using VTKBase: VTKCellData

# function vtk_volume_helper(topo::Topology{D}, 
#     volume::Volume{D},
#     area_vertices::Dict{Int,Vector{Int}}) where D

#     area_ids  = Int[]
#     iterate_volume_areas(area_vertices,topo,volume.id) do area, area_vertices,_
#         push!(area_ids, area.id)
#     end
#     unique_vol_node_ids = get_unique_values(area_vertices,area_ids)
#     return unique_vol_node_ids,get.(Ref(area_vertices),area_ids,nothing)
# end

# """
#     geometry_to_vtk(topo::Topology{D}, filename; verbose=false, save_edges=true) where D

# Write the geometry to a VTK file. In 2D, writes polygonal faces with an `area_id`.
# In 3D, writes polyhedra with a `volume_id`. Additionally, writes a separate
# PolyData file containing all root edges as line cells with an `edge_id` array
# (`<filename>_edges`).
# """
# function geometry_to_vtk(topo::Topology{D}, 
#     filename::String) where {D}


#     raw_points = filter(is_active,get_nodes(topo))

#     points = reduce(hcat,get_coords.(raw_points))

#     if size(points,1) == 2 
#         points = vcat(points,zeros(1,size(points,2)))
#     end


#     node_map = Dict(
#         get_id(node) => i for (i,node) in enumerate(raw_points))



#     edge_cells = [
#         MeshCell(PolyData.Lines(), get.(Ref(node_map),get_edge_node_ids(topo, edge.id),nothing))
#         for edge in RootIterator{D,2}(topo)
#     ]
    


#     area_vertices = Dict(area.id => 
#             get.(Ref(node_map),get_iterative_area_vertex_ids(area, topo),nothing) 
#             for area in RootIterator{D,3}(topo))


#     area_cells = if D ≥ 2 
#         [
#             MeshCell(PolyData.Polys(), area_vertices[area.id])
#             for area in RootIterator{D,3}(topo)
#         ]
#     else
#         nothing
#     end



#     vol_data = Dict(
#         volume.id => vtk_volume_helper(topo,volume,area_vertices)
#         for volume in RootIterator{D,4}(topo)
#     )


#     volume_cells = if D ≥ 3 
#         [
#             VTKPolyhedron(
#                 vol_data[volume.id][1],
#                 vol_data[volume.id][2]...
#             )
#             for volume in RootIterator{D,4}(topo)
#         ]
#     else
#         nothing
#     end






#     saved_files = vtk_multiblock(filename) do vtm  

#         # if u !== nothing

#         #     vtk_grid(vtm,filtered_points ; append = false) do vtk
#         #         vtk["u", VTKPointData()] = u
#         #     end
#         # end

        
#         if D ≥ 3
#             vtk_grid(vtm,points,volume_cells; append = false) do vtk
#                 vtk["volume_id", VTKCellData()] = [volume.id for volume in RootIterator{D,4}(topo)]
#             end
#         end

#         area_block = multiblock_add_block(vtm,filename*"_2d_and_below")
        
#         if D ≥ 2
            
#             vtk_grid(area_block,points,area_cells; append = false) do vtk
#                 vtk["area_id", VTKCellData()] = [area.id for area in RootIterator{D,3}(topo)]
                
#             end
#         end

#         edge_block = multiblock_add_block(area_block,filename*"_edges")


#         vtk_grid(edge_block,points,edge_cells; append = false) do vtk
#             vtk["edge_id", VTKCellData()] = [edge.id for edge in RootIterator{D,2}(topo)]
#         end



#     end
#     nothing
# end

# function geometry_to_vtk(topo::Topology{D}, 
#     filename::String) where {D}


#     raw_points = filter(is_active,get_nodes(topo))

#     points = reduce(hcat,get_coords.(raw_points))

#     if size(points,1) == 2 
#         points = vcat(points,zeros(1,size(points,2)))
#     end

#     # filtered_points = filter(is_active,get_nodes(topo))


#     # active_nodes = filter(is_active,get_nodes(topo))
#     node_map = Dict(
#         get_id(node) => i for (i,node) in enumerate(raw_points))



#     edge_cells = [
#         MeshCell(PolyData.Lines(), get.(Ref(node_map),get_edge_node_ids(topo, edge.id),nothing))
#         for edge in RootIterator{D,2}(topo)
#     ]
    

#     # area_cells = if D ≥ 2 
#     #     [
#     #         MeshCell(PolyData.Polys(), get_area_node_ids(topo, area.id))
#     #         for area in RootIterator{D,3}(topo)
#     #     ]
#     # else
#     #     nothing
#     # end

#     area_vertices = Dict(area.id => 
#             get.(Ref(node_map),get_iterative_area_vertex_ids(area, topo),nothing) 
#             for area in RootIterator{D,3}(topo))


#     area_cells = if D ≥ 2 
#         [
#             MeshCell(PolyData.Polys(), area_vertices[area.id])
#             for area in RootIterator{D,3}(topo)
#         ]
#     else
#         nothing
#     end




#     # volume_cells = if D ≥ 3 
#     #     [
#     #         VTKPolyhedron(
#     #             get_volume_node_ids(topo,volume.id),
#     #             get_area_node_ids(topo)[get_volume_area_ids(topo,volume.id)]...
#     #         ) 
  
#     #     ]
#     # else
#     #     nothing
#     # end
#     vol_data = Dict(
#         volume.id => vtk_volume_helper(topo,volume,area_vertices)
#         for volume in RootIterator{D,4}(topo)
#     )


#     volume_cells = if D ≥ 3 
#         [
#             VTKPolyhedron(
#                 vol_data[volume.id][1],
#                 vol_data[volume.id][2]...
#             )
#             for volume in RootIterator{D,4}(topo)
#         ]
#     else
#         nothing
#     end




#     # extract folder from filename (part until first /)
#     splitted_filename = split(filename,"/")
#     if length(splitted_filename) > 1
#         folder = splitted_filename[1] * "/"
#     else
#         folder = ""
#     end

#     saved_files = vtk_multiblock(filename) do vtm  

#         # if u !== nothing

#         #     vtk_grid(vtm,filtered_points ; append = false) do vtk
#         #         vtk["u", VTKPointData()] = u
#         #     end
#         # end

        
#         if D ≥ 3
#             vtk_grid(vtm,points,volume_cells; append = false) do vtk
#                 vtk["volume_id", VTKCellData()] = [volume.id for volume in RootIterator{D,4}(topo)]
#             end
#         end

#         area_block = multiblock_add_block(vtm,filename*"_2d_and_below")
        
#         if D ≥ 2
            
#             vtk_grid(area_block,points,area_cells; append = false) do vtk
#                 vtk["area_id", VTKCellData()] = [area.id for area in RootIterator{D,3}(topo)]
                
#             end
#         end

#         edge_block = multiblock_add_block(area_block,filename*"_edges")


#         vtk_grid(edge_block,points,edge_cells; append = false) do vtk
#             vtk["edge_id", VTKCellData()] = [edge.id for edge in RootIterator{D,2}(topo)]
#         end



#     end
#     nothing
# end


function vtk_volume_helper(topo::Topology{D}, 
    volume::Volume{D},
    area_vertices) where D

    area_ids  = Int[]
    iterate_volume_areas(area_vertices,topo,volume.id) do area, area_vertices,_
        push!(area_ids, area.id)
    end
    unique_vol_node_ids = get_unique_values(area_vertices,area_ids)
    return unique_vol_node_ids,get.(Ref(area_vertices),area_ids,nothing)
end




@inline trail_zeros3(v::SVector{1,T}) where T= SVector{3}(v[1],zero(T),zero(T))
@inline trail_zeros3(v::SVector{2,T}) where T = SVector{3}(v[1],v[2],zero(T))
@inline trail_zeros3(v::SVector{3,T}) where T = v

function write_vtk(topo::Topology{D}, 
    filename::String = "vtk/sol_to_vtk",
    dh::Union{DofHandler{U},Nothing} = nothing,
    u::Union{AbstractVector,Nothing} = nothing,
    ) where {D,U}

    # either dh and u are provided or not
    @assert (dh !== nothing && u !== nothing) || (dh === nothing && u === nothing)


    raw_points = filter(is_active,get_nodes(topo))



    points = reduce(hcat,get_coords.(raw_points))

    if size(points,1) == 2 
        points = vcat(points,zeros(1,size(points,2)))
    end


    node_map = Dict(
        get_id(node) => i for (i,node) in enumerate(raw_points))
    # node_map = zeros(Int,get_nodes(topo) |> length)
    # for (i,node) in enumerate(raw_points)
    #     node_map[get_id(node)] = i
    # end
    


    edge_cells = [
        MeshCell(PolyData.Lines(), get.(Ref(node_map),get_edge_node_ids(topo, edge.id),nothing))
        # @views MeshCell(PolyData.Lines(), node_map[get_edge_node_ids(topo, edge.id)])
        for edge in RootIterator{D,2}(topo)
    ]
    


    area_vertices = Dict(area.id => 
            get.(Ref(node_map),get_iterative_area_vertex_ids(area, topo),nothing) 
            # @views node_map[get_iterative_area_vertex_ids(area, topo)]
            for area in RootIterator{D,3}(topo))


    area_cells = if D ≥ 2 
        [
            MeshCell(PolyData.Polys(), area_vertices[area.id])
            for area in RootIterator{D,3}(topo)
        ]
    else
        nothing
    end



    vol_data = Dict(
        volume.id => vtk_volume_helper(topo,volume,area_vertices)
        for volume in RootIterator{D,4}(topo)
    )


    volume_cells = if D ≥ 3 
        [
            VTKPolyhedron(
                vol_data[volume.id][1],
                vol_data[volume.id][2]...
            )
            for volume in RootIterator{D,4}(topo)
        ]
    else
        nothing
    end


    if dh !== nothing && u !== nothing
        u_processed = Vector{SVector{3,Float64}}(undef,length(dh.dof_mapping))
        for (n_id,n_dofs) in dh.dof_mapping
            idx = node_map[n_id]

            u_processed[idx] = trail_zeros3(u[n_dofs])
        end
    else
        u_processed = nothing
    end


    vtk_multiblock(filename) do vtm  

        if u_processed !== nothing
            vtk_grid(vtm,points,volume_cells ; append = false) do vtk
                vtk["u", VTKPointData()] = u_processed
            end
        end

        
        if D ≥ 3
            vtk_grid(vtm,points,volume_cells; append = false) do vtk
                vtk["volume_id", VTKCellData()] = [volume.id for volume in RootIterator{D,4}(topo)]
            end
        end

        area_block = multiblock_add_block(vtm,filename*"_2d_and_below")
        
        if D ≥ 2
            
            vtk_grid(area_block,points,area_cells; append = false) do vtk
                vtk["area_id", VTKCellData()] = [area.id for area in RootIterator{D,3}(topo)]
                
            end
        end

        edge_block = multiblock_add_block(area_block,filename*"_edges")


        vtk_grid(edge_block,points,edge_cells; append = false) do vtk
            vtk["edge_id", VTKCellData()] = [edge.id for edge in RootIterator{D,2}(topo)]
        end



    end
    nothing
end