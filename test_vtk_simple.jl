using Ju3VEM

# Test the VTK loader function with a simple test
println("Testing VTK mesh loader...")

try
    # Load the VTK mesh
    mesh = load_vtk_mesh("test_vtk.vtk", StandardEl{1})
    
    println("✓ Successfully loaded VTK mesh")
    println("  - Number of vertices: ", length(get_vertices(mesh)))
    println("  - Number of volumes: ", length(get_volumes(mesh.topo)))
    
    # Print some vertex coordinates
    println("\nFirst 3 vertices:")
    for i in 1:min(3, length(get_vertices(mesh)))
        vertex = get_vertices(mesh)[i]
        println("  Vertex $i: ($(vertex[1]), $(vertex[2]), $(vertex[3]))")
    end
    
    println("\n✓ Test completed successfully!")
    
catch e
    println("✗ Error loading VTK mesh: $e")
    println("Stack trace:")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
end

