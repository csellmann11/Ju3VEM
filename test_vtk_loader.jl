using Ju3VEM

# Test the VTK loader function
println("Testing VTK mesh loader...")

try
    # Load the VTK mesh
    mesh = load_vtk_mesh("test_vtk.vtk", StandardEl{1})
    
    println("✓ Successfully loaded VTK mesh")
    println("  - Number of vertices: ", length(get_vertices(mesh)))
    println("  - Number of volumes: ", length(get_volumes(mesh.topo)))
    println("  - Number of areas: ", length(get_areas(mesh.topo)))
    println("  - Number of edges: ", length(get_edges(mesh.topo)))
    
    # Print some vertex coordinates
    println("\nFirst 5 vertices:")
    for i in 1:min(5, length(get_vertices(mesh)))
        vertex = get_vertices(mesh)[i]
        println("  Vertex $i: ($(vertex[1]), $(vertex[2]), $(vertex[3]))")
    end
    
    # Print volume information
    println("\nVolume information:")
    for (i, volume) in enumerate(RootIterator{4}(mesh.topo))
        println("  Volume $i: ID = $(volume.id)")
    end
    
    println("\n✓ Test completed successfully!")
    
catch e
    println("✗ Error loading VTK mesh: $e")
    rethrow(e)
end

