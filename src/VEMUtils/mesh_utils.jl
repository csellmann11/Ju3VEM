using DelaunayTriangulation: triangulate, voronoi
using Random
# using StableRNGs




function create_rectangular_mesh(
	nx::Int,ny::Int, 
	left::Tuple{Float64, Float64},
	right::Tuple{Float64, Float64},
	::Type{ElT}) where {ElT <: ElType}


	x_coords = range(left[1], stop = right[1], length = nx + 1)
	y_coords = range(left[2], stop = right[2], length = ny + 1)

	coords = [SA[x,y] for x in x_coords, y in y_coords]

	topo = Topology{2}()
	add_node!.(coords, Ref(topo))

	idxs = LinearIndices((nx+1, ny+1))

	for I in CartesianIndices((nx, ny))
		i, j = Tuple(I)
		node_ids = [idxs[i,j],idxs[i+1,j],idxs[i+1,j+1],idxs[i,j+1]]
		add_area!(node_ids,topo)
	end
	return Mesh(topo, ElT())
end

function create_rectangular_mesh(
	nx::Int,ny::Int,
	lx::Float64,
	ly::Float64,
	::Type{ElT}) where {ElT <: ElType}
	
	return create_rectangular_mesh(nx,ny, (0.0,0.0), (lx,ly), ElT)
end

function create_unit_rectangular_mesh(
	nx::Int,ny::Int,
	::Type{ElT}) where {ElT <: ElType}
	
	return create_rectangular_mesh(nx,ny, (0.0,0.0), (1.0,1.0), ElT)
end

function create_rectangular_mesh(
	nx::Int,ny::Int,nz::Int,
	left::Tuple{Float64, Float64, Float64},
	right::Tuple{Float64, Float64, Float64},
	::Type{ElT}) where {ElT <: ElType}
	
	x_coords = range(left[1], stop = right[1], length = nx + 1)
	y_coords = range(left[2], stop = right[2], length = ny + 1)
	z_coords = range(left[3], stop = right[3], length = nz + 1)

	coords = [SA[x,y,z] for x in x_coords, y in y_coords, z in z_coords]
	topo = Topology{3}()
	add_node!.(coords, Ref(topo))

	idxs = LinearIndices((nx+1, ny+1, nz+1))

	for I in CartesianIndices((nx, ny, nz))
		i, j, k = Tuple(I)
		
		# Define face IDs for each hexahedron
		_face_ids = [
			[idxs[i,j,k], idxs[i+1,j,k], idxs[i+1,j+1,k], idxs[i,j+1,k]],           # bottom
			[idxs[i,j,k+1], idxs[i+1,j,k+1], idxs[i+1,j+1,k+1], idxs[i,j+1,k+1]],   # top
			[idxs[i,j,k], idxs[i+1,j,k], idxs[i+1,j,k+1], idxs[i,j,k+1]],           # front 
			[idxs[i+1,j,k], idxs[i+1,j+1,k], idxs[i+1,j+1,k+1], idxs[i+1,j,k+1]],   # right 
			[idxs[i,j,k], idxs[i,j,k+1], idxs[i,j+1,k+1], idxs[i,j+1,k]],           # left
			[idxs[i,j+1,k], idxs[i,j+1,k+1], idxs[i+1,j+1,k+1], idxs[i+1,j+1,k]]    # back
		]
		
		add_volume!(_face_ids, topo)
	end

	return Mesh(topo, ElT())
end

function create_rectangular_mesh(
	nx::Int,ny::Int,nz::Int,
	lx::Float64,ly::Float64,lz::Float64,
	::Type{ElT}) where {ElT <: ElType}
	
	return create_rectangular_mesh(nx,ny,nz, (0.0,0.0), (lx,ly,lz), ElT)
end

function create_unit_rectangular_mesh(
	nx::Int,ny::Int,nz::Int,
	::Type{ElT}) where {ElT <: ElType}
	
	return create_rectangular_mesh(nx,ny,nz, (0.0,0.0,0.0), (1.0,1.0,1.0), ElT)
end

"""
    extrude_to_3d(nz::Int,mesh2d::Mesh{2,ET}) where ET

Extrude a 2D mesh to 3D by adding `nz` layers of nodes in the z-direction. 
The new nodes are added in the order of the edges of the 2D mesh.

# Arguments
- `nz::Int`: The number of layers to add.
- `mesh2d::Mesh{2,ET}`: The 2D mesh to extrude.

# Returns
- `mesh3d::Mesh{3,ET}`: The 3D mesh.
"""
function extrude_to_3d(nz::Int,
    mesh2d::Mesh{2,ET},
    zmax::Float64=1.0) where ET

    topo2d = mesh2d.topo
    vertices = get_vertices(mesh2d)

    topo3    = Topology{3}()

    for vertex in vertices
        z = 0.0 
        add_node!(SA[vertex[1],vertex[2],z], topo3)
    end


    nodes_added = Dict{Int,Int}()

    for element in RootIterator{3}(topo2d)

        old_nodes = get_iterative_area_vertex_ids(element,topo2d)
        new_nodes = Int[]
        for i in 1:nz
            first_node_id = get_area_node_ids(topo2d,element.id)[1]
            face_ids = Vector{Int}[]
            new_z = i * zmax / nz

            local_counter = Ref(1)
            iterate_element_edges(topo2d,element.id) do _, edge_id, _
                # n1id,n2id = get_edge_node_ids(topo2d,edge_id)
                idxp1 = get_next_idx(old_nodes,local_counter[])
                n1id = old_nodes[local_counter[]]
                n2id = old_nodes[idxp1]
                
                n1 = topo3.nodes[n1id]
                n2 = topo3.nodes[n2id]

                n4id = get!(nodes_added,n1id) do  
                    new_coord1 = SA[n1[1],n1[2], new_z]
                    n4id = add_node!(new_coord1, topo3)
                    n4id
                end

                local_counter[] == 1 && push!(new_nodes, n4id)

                
                n3id = get!(nodes_added,n2id) do  
                    new_coord2 = SA[n2[1],n2[2], new_z]
                    n3id = add_node!(new_coord2, topo3)
                    n3id
                end

                idxp1 != 1 && push!(new_nodes, n3id)
                
                push!(face_ids, SA[n1id,n2id,n3id,n4id])
                local_counter[] += 1
            end

            
            push!(face_ids,old_nodes)
            push!(face_ids,new_nodes)
            add_volume!(face_ids, topo3)

            old_nodes .= new_nodes
            empty!(new_nodes)
        end  
    end
    mesh3d = Mesh(topo3, ET())
    return mesh3d
end


function create_voronoi_mesh(left::Tuple{Float64, Float64},
	right::Tuple{Float64, Float64},
	num_points::Int,
	::Type{ElT}, smooth::Bool = true) where {ElT <: ElType}


	len_x = right[1] - left[1]
	len_y = right[2] - left[2]

	nx = floor(Int, len_x / len_y * sqrt(num_points))
	ny = floor(Int, len_y / len_x * nx)

	return create_voronoi_mesh(left, right, nx, ny, ElT, smooth)
end

function create_voronoi_mesh(left::Tuple{Float64, Float64},
	right::Tuple{Float64, Float64},
	nx::Int,ny::Int,
	::Type{ElT}, smooth::Bool = true) where {ElT <: ElType}


	x_coords = range(left[1], stop = right[1], length = nx + 1)
	y_coords = range(left[2], stop = right[2], length = ny + 1)

	L = right[1] - left[1]
	T = right[2] - left[2]


	num_points = nx * ny


	points = Tuple{Float64, Float64}[]

	rng = Random.MersenneTwister(123)
	# add boundary points
	for i in eachindex(x_coords)
		push!(points, (x_coords[i], left[2]))
		push!(points, (x_coords[i], right[2]))
	end
	for i in eachindex(y_coords)
		push!(points, (left[1], y_coords[i]))
		push!(points, (right[1], y_coords[i]))
	end

	n_boundary_points = length(points)

	for i in 1:(num_points)
		push!(points, (rand(rng) * (right[1] - left[1]) + left[1],
			rand(rng) * (right[2] - left[2]) + left[2]))
	end



	unique!(points)
	# Create a Delaunay triangulation
	# rng = StableRNG(123)
	tri = triangulate(points; rng)
	if smooth == false 
		max_iters = 1 
	else 
		max_iters = 10000
	end

	vorn = voronoi(tri, clip = true, smooth = true, rng = rng, maxiters = max_iters)

	topo = Topology{2}()

	id_dict = Dict{Int, Int}()

	for (_, poly_nodes) in vorn.polygons

		poly_nodes = @views poly_nodes[1:end-1]

		poly_nodes = remove_split_edges(poly_nodes, vorn.polygon_points)

		for node_id in poly_nodes
			if haskey(id_dict, node_id)
				continue
			end
			id_dict[node_id] = length(get_nodes(topo)) + 1

			point = vorn.polygon_points[node_id]
			add_node!(SVector(point...), topo)
		end

		id_vec = [id_dict[node_id] for node_id in poly_nodes]

		add_area!(id_vec, topo)
	end

    return Mesh(topo,ElT())
end



"""
	remove_split_edges(poly_nodes::AbstractVector{Int},
	coord_list::AbstractVector{T}) where T

Removes nodes which split a straight edge of the polygon. 
```markdown

5 _____________4 \n
|              | \n
|              | \n
1 ______2______3 \n
```

Here the node #2 splits the edge between #1 and #3. This function removes
node #2 from the polygon.
"""
function remove_split_edges(poly_nodes::AbstractVector{Int},
	coord_list::AbstractVector{T}) where T

	rm = Int[]

	for i in eachindex(poly_nodes)
		ip1 = mod1(i + 1, length(poly_nodes))
		ip2 = mod1(i + 2, length(poly_nodes))

		x1, y1 = coord_list[poly_nodes[i]]
		x2, y2 = coord_list[poly_nodes[ip1]]
		x3, y3 = coord_list[poly_nodes[ip2]]

		area = 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))

		if area < 1e-10
			push!(rm, ip1)
		end
	end

	new_nodes = [poly_nodes[i]
				 for i in eachindex(poly_nodes) if i ∉ rm]

	return new_nodes
end
