using Pkg
Pkg.activate(".")
using Ju3VEM

# Test with a simple VTK file first
println("Testing with simple VTK file...")

try
    mesh = load_vtk_mesh("simple_test.vtk", StandardEl{1})
    println("✓ Simple test successful")
    println("  - Vertices: ", length(get_vertices(mesh)))
    println("  - Volumes: ", length(get_volumes(mesh.topo)))
catch e
    println("✗ Simple test failed: $e")
    rethrow(e)
end
