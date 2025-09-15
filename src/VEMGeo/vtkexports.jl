using WriteVTK 
# using VTKBase: VTKCellData

"""
    geometry_to_vtk(topo::Topology{D}, filename; verbose=false, save_edges=true) where D

Write the geometry to a VTK file. In 2D, writes polygonal faces with an `area_id`.
In 3D, writes polyhedra with a `volume_id`. Additionally, writes a separate
PolyData file containing all root edges as line cells with an `edge_id` array
(`<filename>_edges`).
"""
function geometry_to_vtk(topo::Topology{D}, 
    filename::String,
    u::Union{AbstractVector,Nothing} = nothing) where D


    
    points = reduce(hcat, get_nodes(topo))
    if size(points,1) == 2 
        points = vcat(points,zeros(1,size(points,2)))
    end


    edge_cells = [
        MeshCell(PolyData.Lines(), get_edge_node_ids(topo, edge.id))
        for edge in RootIterator{D,2}(topo)
    ]
    

    area_cells = if D ≥ 2 
        [
            MeshCell(PolyData.Polys(), get_area_node_ids(topo, area.id))
            for area in RootIterator{D,3}(topo)
        ]
    else
        nothing
    end



    volume_cells = if D ≥ 3 
        [
            VTKPolyhedron(
                get_volume_node_ids(topo,volume.id),
                get_area_node_ids(topo)[get_volume_area_ids(topo,volume.id)]...
            ) 
            for volume in RootIterator{D,4}(topo)
        ]
    else
        nothing
    end


    # uSA = [SA[u[i],0.0,0.0] for i in eachindex(u)]
    # point_coords = get_coords.(topo.nodes)



    # extract folder from filename (part until first /)
    splitted_filename = split(filename,"/")
    if length(splitted_filename) > 1
        folder = splitted_filename[1] * "/"
    else
        folder = ""
    end
    



    saved_files = vtk_multiblock(filename) do vtm  

        if u !== nothing

            vtk_grid(vtm,points,volume_cells; append = false) do vtk
                vtk["u", VTKPointData()] = u
            end
        end

        
        if D ≥ 3
            vtk_grid(vtm,points,volume_cells; append = false) do vtk
                vtk["volume_id", VTKCellData()] = [volume.id for volume in RootIterator{D,4}(topo)]
            end
        end

        area_block = multiblock_add_block(vtm,folder*"2d_and_below")
        
        if D ≥ 2
            
            vtk_grid(area_block,points,area_cells; append = false) do vtk
                vtk["area_id", VTKCellData()] = [area.id for area in RootIterator{D,3}(topo)]
                
            end
        end

        edge_block = multiblock_add_block(area_block,folder*"edges")


        vtk_grid(edge_block,points,edge_cells; append = false) do vtk
            vtk["edge_id", VTKCellData()] = [edge.id for edge in RootIterator{D,2}(topo)]
        end



    end
    nothing
end