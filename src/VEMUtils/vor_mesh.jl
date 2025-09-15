using DelaunayTriangulation: triangulate, voronoi
using Random
# using StableRNGs

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
