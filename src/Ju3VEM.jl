module Ju3VEM

using Reexport

module VEMGeo
    # Core dependencies commonly used across files
    using StaticArrays
    using LinearAlgebra
    using Random
    using Statistics
    using FixedSizeArrays
    using OrderedCollections
    using SmallCollections
    using WriteVTK
    using LoopVectorization
    using Collects
    using Symbolics
    using Bumper
    using FastGaussQuadrature
    using LazyArrays
    # Symbolics.RuntimeGeneratedFunctions.init(@__MODULE__)

    # Load core geometry/topology and utilities first (dependency order matters)
    include("VEMGeo/topo.jl")
    include("VEMGeo/small_utils.jl")
    include("VEMGeo/lagrange_utils.jl")

    # Refinement / coarsening and I/O
    include("VEMGeo/element_refinement.jl")
    include("VEMGeo/element_coarsening.jl")
    include("VEMGeo/vtkexports.jl")

    # Polynomials (includes monomials and stretched_matrices internally)
    include("VEMGeo/Polynomials/polynomials.jl")
    include("VEMGeo/Polynomials/poly_tansform.jl")
    include("VEMGeo/sym_integration.jl")

    # Triangulation and 3D utilities (flattened; no submodules inside the file)
    include("VEMGeo/triangulation.jl")
    include("VEMGeo/integration.jl")
    include("flatten_vecs.jl")

    # include("mesh.jl")
    # include("face_projector.jl")
    # include("volume_projector.jl")

    # Exports: expose primary types and functions
    export 
        # Types
        Node, NManifold, Topology, Edge, Area, Volume, RootIterator, 
        Polynomial, Monomial, PolynomialBase, BaseInfo, StretchedMatrix, 
        TriangleQuadRule, Topology, FlattenVecs,
        # Topology getters/setters
        get_coords, get_id, is_active, is_root, 
        get_nodes, get_edges, get_areas, get_volumes, 
        get_volume_node_ids, get_volume_edge_ids, get_volume_area_ids, 
        get_area_node_ids, get_area_edge_ids, get_edge_node_ids, 
        add_node!, add_edge!, add_area!, add_volume!, 
        get_iterative_area_vertex_ids, iterate_element_edges, iterate_volume_areas, 
        apply_f_on, apply_f_on_roots, num_roots, 
        # Utils
        get_next_idx, get_prev_idx, get_unique_values, find_single_intersec, max_node_distance,
        # Refinement / Coarsening
        _refine!, _coarsen!, downstream_deactivate!, 
        # VTK Export
        geometry_to_vtk, 
        # Lagrange utilities
        lagrange_shape_function, interpolate_edge, 
        # Monomial API
        derivative, ∇, ∂, 
        # Polynomials API
        get_base, get_bi_base, sol_proj, grad_sol_proj, get_ϕi, ∇p, ∇x, poly_grad, 
        # Triangulation / 3D helpers
        triangulate_polygon, triangulate_planar_polygon3D, 
        tetrahedralize_points_convex, tetrahedralize_volume_local_ids, 
        tetrahedralize_volume, build_tet_topology_from_volume, 
        FaceTriangulations3D, get_area_triangles, polygon_area3D, 
        integrate_polynomial_over_face, integrate_polynomial_over_volume, 
        volume_from_faces, volume_of_topo_volume,
        # Poly transform
        compute_transformation_coeffs2d_to_3d, compute_transformation_coeffs3d_to_2d,
        compute_transformation_coeffs3d_to_2d!,compute_transformation_coeffs2d_to_2d!,
        compute_transformation_coeffs2d_to_2d,
        # Sym integration
        precompute_face_monomials, compute_face_integral, compute_face_integral_unshifted,
        get_area, get_bc
end # module VEMGeo

# Re-export the VEMGeo API from the top-level module
@reexport using .VEMGeo

end # module Ju3VEM
