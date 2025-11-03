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
	
	return create_rectangular_mesh(nx,ny,nz, (0.0,0.0,0.0), (lx,ly,lz), ElT)
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


	x_coords = range(left[1], stop = right[1], length = 2)
	y_coords = range(left[2], stop = right[2], length = 2)

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

		if area < 1e-13
			push!(rm, ip1)
		end
	end

	new_nodes = [poly_nodes[i]
				 for i in eachindex(poly_nodes) if i ∉ rm]

	return new_nodes
end

"""
    load_vtk_mesh(filename::String, ::Type{ElT}) where {ElT <: ElType}

Load a VTK unstructured grid file and convert it to a custom Mesh structure.

# Arguments
- `filename::String`: Path to the VTK file
- `::Type{ElT}`: Element type (e.g., StandardEl{1}, SerendipityEl{1})

# Returns
- `mesh::Mesh{3,ElT}`: 3D mesh structure

# Example
```julia
mesh = load_vtk_mesh("mesh.vtk", StandardEl{1})
```
"""
function load_vtk_mesh(filename::String, ::Type{ElT}) where {ElT <: ElType}
    # Read the VTK file
    lines = readlines(filename)
    
    # Find the POINTS section
    points_start = findfirst(line -> occursin("POINTS", line), lines)
    @assert points_start !== nothing "POINTS section not found in VTK file"
    
    # Parse number of points
    points_line = lines[points_start]
    n_points = parse(Int, split(points_line)[2])
    
    # Read point coordinates
    points = Vector{SVector{3,Float64}}()
    for i in (points_start+1):(points_start+n_points)
        coords = parse.(Float64, split(lines[i]))
        push!(points, SVector{3,Float64}(coords[1], coords[2], coords[3]))
    end
    
    # Find the CELLS section
    cells_start = findfirst(line -> occursin("CELLS", line), lines)
    @assert cells_start !== nothing "CELLS section not found in VTK file"
    
    # Parse number of cells and total connectivity entries
    cells_line = lines[cells_start]
    n_cells, total_connectivity = parse.(Int, split(cells_line)[2:3])
    
    # Read cell connectivity
    cells = Vector{Vector{Vector{Int}}}()  # Each cell contains multiple faces
    line_idx = cells_start + 1
    
    for i in 1:n_cells
        # Read all data for this cell (may span multiple lines)
        cell_data = Int[]
        current_line = lines[line_idx]
        line_data = parse.(Int, split(current_line))
        total_count = line_data[1]  # Total count of data entries (ignore for now)
        n_faces = line_data[2]      # Number of faces
        
        # Add the rest of the first line (skip total_count and n_faces)
        append!(cell_data, line_data[3:end])
        
        # Read additional lines if needed - continue until we have enough data
        while true
            # Calculate how much data we need: n_faces + sum of all face vertex counts
            needed_data = n_faces
            data_idx = 1
            for face_idx in 1:n_faces
                if data_idx > length(cell_data)
                    break
                end
                n_face_vertices = cell_data[data_idx]
                needed_data += n_face_vertices
                data_idx += n_face_vertices + 1
            end
            
            if length(cell_data) >= needed_data
                break
            end
            
            # Need more data - read next line
            line_idx += 1
            if line_idx > length(lines)
                error("Unexpected end of file while parsing cell $i")
            end
            current_line = lines[line_idx]
            line_data = parse.(Int, split(current_line))
            append!(cell_data, line_data)
        end
        
        
        # Parse faces for this cell
        cell_faces = Vector{Vector{Int}}()
        data_idx = 1
        
        for face_idx in 1:n_faces
            if data_idx > length(cell_data)
                error("Not enough data for face $face_idx of cell $i")
            end
            
            n_face_vertices = cell_data[data_idx]
            if data_idx + n_face_vertices > length(cell_data)
                error("Face $face_idx data incomplete for cell $i")
            end
            
            face_vertices = cell_data[(data_idx+1):(data_idx+n_face_vertices)]
            push!(cell_faces, face_vertices)
            data_idx += n_face_vertices + 1
        end
        
        push!(cells, cell_faces)
        line_idx += 1
    end
    
    # Find the CELL_TYPES section
    types_start = findfirst(line -> occursin("CELL_TYPES", line), lines)
    @assert types_start !== nothing "CELL_TYPES section not found in VTK file"
    
    # Read cell types
    cell_types = Int[]
    for i in (types_start+1):(types_start+n_cells)
        push!(cell_types, parse(Int, lines[i]))
    end
    
    # Create topology
    topo = Topology{3}()
    
    # Add nodes to topology
    for (i, point) in enumerate(points)
        add_node!(point, topo, i)
    end
    
    # Process cells (assuming all are polyhedra - VTK cell type 42)
    for (cell_idx, cell_faces) in enumerate(cells)
        if cell_types[cell_idx] == 42  # VTK_POLYHEDRON
            # Convert face vertex indices to 1-based indexing (VTK uses 0-based)
            faces_1based = [face_vertices .+ 1 for face_vertices in cell_faces]
            
            # Add volume with faces
            add_volume!(faces_1based, topo)
        else
            @warn "Unsupported cell type $(cell_types[cell_idx]) for cell $cell_idx"
        end
    end
    
    # Create mesh
    return Mesh(topo, ElT())
end



# ===================== 3D Voronoi (brick-bounded, no external libs) =====================

@inline _as_SA(p::Tuple{Float64,Float64,Float64}) = SA[p[1],p[2],p[3]]
@inline _as_SA(p::SVector{3,Float64}) = p

@inline function _plane_eval(n::SVector{3,Float64}, d::Float64, x::SVector{3,Float64})
    return dot(n, x) - d
end

@inline function _line_plane_intersection(p1::SVector{3,Float64}, p2::SVector{3,Float64}, n::SVector{3,Float64}, d::Float64)
    denom = dot(n, p2 - p1)
    t = (d - dot(n, p1)) / denom
    return p1 + t * (p2 - p1)
end

@inline function _orthonormal_basis(n::SVector{3,Float64})
    a = abs(n[1]) < 0.9 ? SA[1.0,0.0,0.0] : SA[0.0,1.0,0.0]
    u = normalize(cross(n, a))
    v = cross(n, u)
    return u, v
end

function _order_points_on_plane(points::Vector{SVector{3,Float64}}, n::SVector{3,Float64})
    c = reduce(+, points) / length(points)
    u, v = _orthonormal_basis(n)
    angles = map(p -> atan(dot(p - c, v), dot(p - c, u)), points)
    idx = sortperm(angles)
    return points[idx]
end

@inline function _is_close(a::SVector{3,Float64}, b::SVector{3,Float64}, tol::Float64)
    return maximum(abs.(a .- b)) <= tol
end

function _dedup_points(points::Vector{SVector{3,Float64}}, tol::Float64)
    # Deduplicate with tolerance using quantization key
    q = (p)->(round(Int, p[1]/tol), round(Int, p[2]/tol), round(Int, p[3]/tol))
    seen = Dict{NTuple{3,Int},Int}()
    out = SVector{3,Float64}[]
    for p in points
        key = q(p)
        if haskey(seen, key)
            # Optionally refine by exact tol check against existing
            existing = out[seen[key]]
            if !_is_close(existing, p, tol)
                # Rare quantization clash; keep distinct by slight key perturbation
                push!(out, p)
            end
        else
            seen[key] = length(out) + 1
            push!(out, p)
        end
    end
    return out
end

function _dedup_consecutive(poly::Vector{SVector{3,Float64}}, tol::Float64)
    isempty(poly) && return poly
    out = SVector{3,Float64}[]
    last = poly[1]
    push!(out, last)
    for i in 2:length(poly)
        p = poly[i]
        if !_is_close(p, last, tol)
            push!(out, p)
            last = p
        end
    end
    # also check first/last
    if length(out) >= 2 && _is_close(out[1], out[end], tol)
        pop!(out)
    end
    return out
end

function _clip_polygon_by_plane(poly::Vector{SVector{3,Float64}}, n::SVector{3,Float64}, d::Float64, tol::Float64)
    # Sutherland–Hodgman polygon clip in 3D against plane dot(n,x) <= d
    newpoly = SVector{3,Float64}[]
    added = SVector{3,Float64}[]  # intersections to form cap face later
    if isempty(poly)
        return newpoly, added
    end
    s_prev = _plane_eval(n, d, poly[end])
    prev = poly[end]
    for curr in poly
        s_curr = _plane_eval(n, d, curr)
        inside_prev = s_prev <= tol
        inside_curr = s_curr <= tol
        if inside_curr
            if inside_prev
                # in->in
                push!(newpoly, curr)
            else
                # out->in
                inter = _line_plane_intersection(prev, curr, n, d)
                push!(newpoly, inter)
                push!(newpoly, curr)
                push!(added, inter)
            end
        else
            if inside_prev
                # in->out
                inter = _line_plane_intersection(prev, curr, n, d)
                push!(newpoly, inter)
                push!(added, inter)
            else
                # out->out -> nothing
            end
        end
        prev = curr
        s_prev = s_curr
    end
    newpoly = _dedup_consecutive(newpoly, tol)
    return newpoly, added
end

function _clip_polyhedron_by_plane(faces::Vector{Vector{SVector{3,Float64}}}, n::SVector{3,Float64}, d::Float64, tol::Float64)
    new_faces = Vector{Vector{SVector{3,Float64}}}()
    ring = SVector{3,Float64}[]
    for face in faces
        clipped, added = _clip_polygon_by_plane(face, n, d, tol)
        if length(clipped) >= 3
            push!(new_faces, clipped)
        end
        append!(ring, added)
    end
    # Add cap face if the plane intersected the polyhedron
    if !isempty(ring)
        ring = _dedup_points(ring, max(tol, 10*eps(Float64)))
        if length(ring) >= 3
            ring = _order_points_on_plane(ring, n)
            push!(new_faces, ring)
        end
    end
    return new_faces
end

function _brick_faces(left::SVector{3,Float64}, right::SVector{3,Float64})
    x0,y0,z0 = left
    x1,y1,z1 = right
    c000 = SA[x0,y0,z0]; c100 = SA[x1,y0,z0]; c110 = SA[x1,y1,z0]; c010 = SA[x0,y1,z0]
    c001 = SA[x0,y0,z1]; c101 = SA[x1,y0,z1]; c111 = SA[x1,y1,z1]; c011 = SA[x0,y1,z1]
    return Vector{SVector{3,Float64}}[
        [c000,c100,c110,c010], # z = z0 (bottom)
        [c001,c011,c111,c101], # z = z1 (top)
        [c000,c010,c011,c001], # x = x0 (left)
        [c100,c101,c111,c110], # x = x1 (right)
        [c000,c001,c101,c100], # y = y0 (front)
        [c010,c110,c111,c011], # y = y1 (back)
    ]
end

@inline function _bisector_plane(si::SVector{3,Float64}, sj::SVector{3,Float64})
    n = sj - si
    d = (dot(sj, sj) - dot(si, si)) / 2
    return n, d
end

function _create_cell_faces_for_seed(seed::SVector{3,Float64}, seeds::AbstractVector{SVector{3,Float64}}, left::SVector{3,Float64}, right::SVector{3,Float64}, tol::Float64)
    faces = _brick_faces(left, right)
    for s in seeds
        if s === seed; continue; end
        n, d = _bisector_plane(seed, s)
        faces = _clip_polyhedron_by_plane(faces, n, d, tol)
        if isempty(faces)
            break
        end
    end
    return faces
end

function _add_polyhedron_to_topology!(faces::Vector{Vector{SVector{3,Float64}}}, topo::Topology{3}, dedup_tol::Float64)
    # Global dedup by quantization map
    quant = (p)->(round(Int, p[1]/dedup_tol), round(Int, p[2]/dedup_tol), round(Int, p[3]/dedup_tol))
    # Reuse existing nodes already present in topo
    node_map = Dict{NTuple{3,Int},Int}()
    # seed with current nodes in topology
    for n in get_nodes(topo)
        key = quant(SA[n[1],n[2],n[3]])
        node_map[key] = get_id(n)
    end
    # build faces in terms of node ids
    face_ids_col = Vector{Vector{Int}}()
    for poly in faces
        poly = _dedup_consecutive(poly, dedup_tol)
        if length(poly) < 3; continue; end
        ids = Int[]
        for p in poly
            key = quant(p)
            id = get!(node_map, key) do
                add_node!(p, topo)
            end
            push!(ids, id)
        end
        # remove duplicate last==first if present
        if !isempty(ids) && first(ids) == last(ids)
            pop!(ids)
        end
        if length(ids) >= 3
            push!(face_ids_col, ids)
        end
    end
    if !isempty(face_ids_col)
        add_volume!(face_ids_col, topo)
    end
end

"""
	create_voronoi_mesh_3d(left, right, seeds, ::Type{ElT}; dedup_tol=defaultTol)

Generate a 3D Voronoi tessellation clipped to the axis-aligned brick defined by
`left = (x0,y0,z0)` and `right = (x1,y1,z1)`, using the provided seed points.

- No external meshing libraries are used.
- Vertices are deduplicated with tolerance `dedup_tol` to avoid near-duplicates.
- Returns a `Mesh{3,ElT}` ready for VTK export via `write_vtk`.
"""
function create_voronoi_mesh_3d(
	left::Tuple{Float64,Float64,Float64},
	right::Tuple{Float64,Float64,Float64},
	seeds_in::AbstractVector{<:SVector{3,Float64}},
	::Type{ElT}; dedup_tol::Float64 = -1.0,
) where {ElT<:ElType}

	left3 = _as_SA(left)
	right3 = _as_SA(right)
	seeds = collect(seeds_in)

	# Default dedup tolerance based on domain diagonal
	default_tol = norm(right3 - left3) * 1e-10
	tol = default_tol
	ddtol = dedup_tol < 0 ? default_tol : dedup_tol

	topo = Topology{3}()

	for i in eachindex(seeds)
		faces = _create_cell_faces_for_seed(seeds[i], seeds, left3, right3, tol)
		isempty(faces) && continue
		_add_polyhedron_to_topology!(faces, topo, ddtol)
	end

	return Mesh(topo, ElT())
end

"""
	create_voronoi_mesh_3d(left, right, num_points, ::Type{ElT}; rng=MersenneTwister(123), dedup_tol=defaultTol)

Convenience wrapper that samples `num_points` uniform seeds in the brick and
builds the Voronoi mesh.
"""
function create_voronoi_mesh_3d(
	left::Tuple{Float64,Float64,Float64},
	right::Tuple{Float64,Float64,Float64},
	num_points::Int,
	::Type{ElT}; rng = Random.MersenneTwister(123), dedup_tol::Float64 = -1.0,
) where {ElT<:ElType}

	left3 = _as_SA(left)
	right3 = _as_SA(right)
	L = right3 - left3
	seeds = [SA[left3[1] + rand(rng)*L[1], left3[2] + rand(rng)*L[2], left3[3] + rand(rng)*L[3]] for _ in 1:num_points]
	return create_voronoi_mesh_3d(left, right, seeds, ElT; dedup_tol)
end

"""
	write_voronoi3d_vtk(left, right, seeds_or_n, filename, ::Type{ElT}=StandardEl{1})

Helper to generate and immediately write a Voronoi 3D mesh to VTK using the
existing `write_vtk` utilities.
"""
function write_voronoi3d_vtk(
	left::Tuple{Float64,Float64,Float64},
	right::Tuple{Float64,Float64,Float64},
	seeds::AbstractVector{<:SVector{3,Float64}},
	filename::String,
	::Type{ElT}=StandardEl{1}; dedup_tol::Float64 = -1.0,
) where {ElT<:ElType}

	mesh = create_voronoi_mesh_3d(left, right, seeds, ElT; dedup_tol)
	write_vtk(mesh.topo, filename)
	return mesh
end

function write_voronoi3d_vtk(
	left::Tuple{Float64,Float64,Float64},
	right::Tuple{Float64,Float64,Float64},
	num_points::Int,
	filename::String,
	::Type{ElT}=StandardEl{1}; rng = Random.MersenneTwister(123), dedup_tol::Float64 = -1.0,
) where {ElT<:ElType}

	mesh = create_voronoi_mesh_3d(left, right, num_points, ElT; rng, dedup_tol)
	write_vtk(mesh.topo, filename)
	return mesh
end


# ===================== Lloyd relaxation for 3D Voronoi =====================

@inline function _tet_volume(a::SVector{3,Float64}, b::SVector{3,Float64}, c::SVector{3,Float64}, d::SVector{3,Float64})
    return abs(dot(cross(b - a, c - a), d - a)) / 6
end

function _convex_poly_centroid(points::Vector{SVector{3,Float64}})
    # Triangulate convex hull of points and compute volume centroid
    tets = tetrahedralize_points_convex(points)
    Vsum = 0.0
    Csum = SA[0.0, 0.0, 0.0]
    for (i,j,k,l) in tets
        a = points[i]; b = points[j]; c = points[k]; d = points[l]
        V = _tet_volume(a,b,c,d)
        Vsum += V
        Csum += V * ((a + b + c + d) / 4)
    end
    if Vsum == 0.0
        # Fallback: mean of vertices
        return reduce(+, points) / length(points)
    end
    return Csum / Vsum
end

function _cell_centroid_from_faces(faces::Vector{Vector{SVector{3,Float64}}})
    # Gather unique vertices
    pts = SVector{3,Float64}[]
    seen = Set{NTuple{3,Float64}}()
    for f in faces
        for p in f
            t = (p[1],p[2],p[3])
            if !(t in seen)
                push!(seen, t)
                push!(pts, p)
            end
        end
    end
    return _convex_poly_centroid(pts)
end

"""
	relax_voronoi3d_seeds(left, right, seeds; maxiters=10, move_tol=1e-6, step=1.0, min_spacing=nothing, spacing_passes=1, spacing_step=1.0)

Lloyd-style relaxation for 3D Voronoi seeds with optional spacing control.

- `step` (0,1]: under-relaxation for centroid move.
- `min_spacing`: if set (e.g., 0.02*min(box size)), applies a min-distance projection after each Lloyd step to reduce very short edges.
- `spacing_passes`: number of projection sweeps per iteration (1–2 is usually enough).
- `spacing_step` (0,1]: fraction of the computed projection applied per pass (for stability).
"""
function relax_voronoi3d_seeds(
	left::Tuple{Float64,Float64,Float64},
	right::Tuple{Float64,Float64,Float64},
	seeds_in::AbstractVector{<:SVector{3,Float64}};
	maxiters::Int=10, move_tol::Float64=1e-6, step::Float64=1.0,
	min_spacing::Union{Nothing,Float64}=nothing, spacing_passes::Int=1, spacing_step::Float64=1.0,
    verbose::Bool=false
)::Vector{SVector{3,Float64}}

	left3 = _as_SA(left)
	right3 = _as_SA(right)
	L = right3 - left3
	seeds = collect(seeds_in)

	# tolerance for clipping geometry (robustness)
	clip_tol = norm(L) * 1e-10

	for it in 1:maxiters
		new_seeds = similar(seeds)
		dmax = 0.0
		for i in eachindex(seeds)
			faces = _create_cell_faces_for_seed(seeds[i], seeds, left3, right3, clip_tol)
			if isempty(faces)
				new_seeds[i] = seeds[i] # keep if degenerate
				continue
			end
			c = _cell_centroid_from_faces(faces)
			# bounded step and clamp to domain
			ns = (1 - step) * seeds[i] + step * c
			ns = SA[
				clamp(ns[1], left3[1], right3[1]),
				clamp(ns[2], left3[2], right3[2]),
				clamp(ns[3], left3[3], right3[3])
			]
			new_seeds[i] = ns
			dmax = max(dmax, maximum(abs.(ns .- seeds[i])) / maximum(L))
		end

		# Optional min-distance projection to improve shortest/longest edge ratio
		if min_spacing !== nothing && min_spacing > 0
			ms2 = min_spacing^2
			for _p in 1:max(1, spacing_passes)
				delta = [SA[0.0,0.0,0.0] for _ in eachindex(new_seeds)]
				for i in 1:length(new_seeds)-1
					pi = new_seeds[i]
					for j in (i+1):length(new_seeds)
						pj = new_seeds[j]
						dv = pj - pi
						d2 = sum(abs2, dv)
						if d2 < ms2
							if d2 == 0.0
								# split deterministically along x
								dir = SA[1.0,0.0,0.0]
								gap = min_spacing
								shift = 0.5 * spacing_step * gap * dir
								delta[i] -= shift
								delta[j] += shift
							else
								d = sqrt(d2)
								gap = min_spacing - d
								dir = dv / d
								shift = 0.5 * spacing_step * gap * dir
								delta[i] -= shift
								delta[j] += shift
							end
						end
					end
				end
				# Apply accumulated shifts and clamp to domain
				for i in eachindex(new_seeds)
					p = new_seeds[i] + delta[i]
					new_seeds[i] = SA[
						clamp(p[1], left3[1], right3[1]),
						clamp(p[2], left3[2], right3[2]),
						clamp(p[3], left3[3], right3[3])
					]
				end
			end
		end
		seeds = new_seeds
        if verbose
            println("Iteration $it, max move: $dmax")
        end
		if dmax < move_tol
			break
		end
	end

	return seeds
end

function relax_voronoi3d(
	left::Tuple{Float64,Float64,Float64},
	right::Tuple{Float64,Float64,Float64},
	seeds::AbstractVector{<:SVector{3,Float64}},
	::Type{ElT}; maxiters::Int=10, move_tol::Float64=1e-6, step::Float64=1.0,
	min_spacing::Union{Nothing,Float64}=nothing, spacing_passes::Int=1, spacing_step::Float64=1.0, dedup_tol::Float64=-1.0,
) where {ElT<:ElType}

	seeds_relaxed = relax_voronoi3d_seeds(left, right, seeds; maxiters, move_tol, step, min_spacing, spacing_passes, spacing_step)
	return create_voronoi_mesh_3d(left, right, seeds_relaxed, ElT; dedup_tol)
end