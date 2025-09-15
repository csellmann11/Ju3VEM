using Ju3VEM




coords = SVector{3, Float64}[
    [0.6875, 0.0, 0.0], 
    [0.75, 0.0, 0.0], 
    [0.6875, 0.0625, 0.0],
    [0.75, 0.0625, 0.0],  

    [0.6875, 0.0, 0.0625], 
    [0.75, 0.0, 0.0625], 
    [0.6796875, 0.0546875, 0.0703125],
    [0.7578125, 0.0703125, 0.0703125], 
    ]


coords = SVector{3, Float64}[
    [0.6875, 0.0, 0.0], 
    [0.75, 0.0, 0.0], 
    [0.6875, 0.0625, 0.0], 
    [0.75, 0.0625, 0.0], 
  

    [0.6875, 0.0, 0.0625], 
    [0.75, 0.0, 0.0625], 
    [0.6875, 0.0625, 0.0625],
    [0.75+0.04, 0.0625, 0.0625], 
]


# Initialize topology
topo = Topology{3}()

# Add all nodes to topology
add_node!.(coords, Ref(topo))
idxs = LinearIndices((2,2,2))
for I in CartesianIndices((1,1,1))
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

# =======
mesh = Mesh(topo, StandardEl{1}())
cv = CellValues{1}(mesh)

element1 = get_volumes(topo)[1] 

reinit!(element1.id,cv)

cv.volume_data.integrals
# get_field.(cv.facedata_col, :dΩ)

[cv.facedata_col[i].dΩ.integrals[1] for i in 1:length(cv.facedata_col)]

get_volume_area_ids(topo, element1.id)


function face_area(vertices)
    n = length(vertices)
    area = 0.0
    for i in 2:n-1
        v1 = vertices[i] - vertices[1]
        v2 = vertices[i+1] - vertices[1]
        area += norm(cross(v1, v2)) / 2
    end
    return area
end


for area_id in get_volume_area_ids(topo, element1.id)
    node_ids = get_area_node_ids(topo, area_id)
    vertices = get_nodes(topo)[node_ids]
    area_id == 4 && @show vertices
    area = face_area(vertices)
    @show area_id, area

end