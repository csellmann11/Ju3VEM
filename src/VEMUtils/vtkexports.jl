
function vtk_volume_helper(topo::Topology{D}, 
    volume::Volume{D},
    area_vertices) where D

    area_ids  = Int32[]
    iterate_volume_areas(area_vertices,topo,volume.id) do area, area_vertices,_
        push!(area_ids, area.id)
    end
    # unique_vol_node_ids = get_unique_values(area_vertices,area_ids)
    vertex_ids = get.(Ref(area_vertices),area_ids,nothing)
    unique_vol_node_ids = unique(Iterators.flatten(vertex_ids))
    return unique_vol_node_ids,vertex_ids
end




@inline trail_zeros3(v::StaticVector{1,T}) where T = SVector{3}(v[1],zero(T),zero(T))
@inline trail_zeros3(v::StaticVector{2,T}) where T = SVector{3}(v[1],v[2],zero(T))
@inline trail_zeros3(v::StaticVector{3,T}) where T = SVector{3}(v)

function write_vtm(topo::Topology{D}, 
    filename::String = "vtk/sol_to_vtk",
    dh::Union{DofHandler{U},Nothing} = nothing,
    u::Union{AbstractVector,Nothing} = nothing;
    cell_data = nothing
    ) where {D,U}

    # either dh and u are provided or not
    @assert (dh !== nothing && u !== nothing) || (dh === nothing && u === nothing)


    raw_points = filter(is_active,get_nodes(topo))



    points = reduce(hcat,get_coords.(raw_points))

    if size(points,1) == 2 
        points = vcat(points,zeros(1,size(points,2)))
    end



    node_map = zeros(Int,get_nodes(topo) |> length)
    for (i,node) in enumerate(raw_points)
        node_map[get_id(node)] = i
    end
    


    edge_cells = [
        @views MeshCell(PolyData.Lines(), node_map[get_edge_node_ids(topo, edge.id)])
        for edge in RootIterator{2}(topo)
    ]
    


    area_vertices = Dict(area.id => 
            @views node_map[get_iterative_area_vertex_ids(area, topo)]
            for area in RootIterator{3}(topo))




    area_cells = if D ≥ 2 
        [
            MeshCell(PolyData.Polys(), area_vertices[area.id])
            for area in RootIterator{3}(topo)
        ]
    else
        nothing
    end



    vol_data = Dict(
        volume.id => vtk_volume_helper(topo,volume,area_vertices)
        for volume in RootIterator{4}(topo)
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

        if cell_data !== nothing
            vtk_grid(vtm,points,volume_cells; append = false) do vtk
                vtk["cell_data", VTKCellData()] = cell_data
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



function write_vtk(topo::Topology{D}, 
    filename::String = "vtk/sol_to_vtk",
    dh::Union{DofHandler{U},Nothing} = nothing,
    u::Union{AbstractVector,Nothing} = nothing;
    cell_data_col::Tuple = ()
    ) where {D,U}

    # either dh and u are provided or not
    @assert (dh !== nothing && u !== nothing) || (dh === nothing && u === nothing)


    raw_points = filter(x -> is_active(x,topo),get_nodes(topo))



    points = reduce(hcat,get_coords.(raw_points))

    if size(points,1) == 2 
        points = vcat(points,zeros(1,size(points,2)))
    end



    node_map = zeros(Int,get_nodes(topo) |> length)
    for (i,node) in enumerate(raw_points)
        node_map[get_id(node)] = i
    end
    




    area_vertices = Dict(area.id => 
            @views node_map[get_iterative_area_vertex_ids(area, topo)]
            for area in RootIterator{3}(topo))


    vol_data = Dict(
        volume.id => vtk_volume_helper(topo,volume,area_vertices)
        for volume in RootIterator{4}(topo)
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

    vtk_grid(filename,points,volume_cells ; append = false) do vtk
        if u_processed !== nothing
        vtk["u", VTKPointData()] = u_processed
        end
        
        if !isempty(cell_data_col)
            for i in eachindex(cell_data_col)
                vtk["cell_data_$i", VTKCellData()] = cell_data_col[i]
            end
        end

        if D == 3 
            vtk["volume_id", VTKCellData()] = [volume.id for volume in RootIterator{D,4}(topo)]
        elseif D == 2
            vtk["area_id", VTKCellData()] = [area.id for area in RootIterator{D,3}(topo)]
        else
            throw(ErrorException("Unsupported dimension: $D"))
        end
    end
    
    nothing
end