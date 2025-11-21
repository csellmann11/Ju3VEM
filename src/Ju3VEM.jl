module Ju3VEM

using Reexport

include("StructMechanics/material_laws.jl")
@reexport using .MaterialLaws


import Tensorial 
using Tensorial: ⊡₂
using StaticArrays

function Tensorial.:⊡₂(a::V1,b::V2) where {V1<:StaticArray,V2<:StaticArray}
    res = Tensorial.Tensor(a) ⊡₂ Tensorial.Tensor(b)

    if res isa Tensorial.Tensor
        return SArray{Tuple{size(res)...}}(res.data)
    else
        return res
    end
end

export ⊡₂

module VEMGeo
    # Core dependencies commonly used across files
    using StaticArrays
    using LinearAlgebra
    using Octavian
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
    using SparseArrays

    # Symbolics.RuntimeGeneratedFunctions.init(@__MODULE__)

    include("flatten_vecs.jl")
    # Load core geometry/topology and utilities first (dependency order matters)
    include("VEMGeo/topo.jl")
    include("VEMGeo/small_utils.jl")
    include("VEMGeo/lagrange_utils.jl")

    # Refinement / coarsening and I/O
    include("VEMGeo/element_refinement.jl")
    include("VEMGeo/element_coarsening.jl")
    # include("VEMGeo/vtkexports.jl")


    # include("element_mapping.jl")#! relies on eltye and mesh

    # Polynomials (includes monomials and stretched_matrices internally)
    include("VEMGeo/Polynomials/polynomials.jl")
    include("VEMGeo/Polynomials/poly_tansform.jl")
    include("VEMGeo/sym_integration.jl")

    # Triangulation and 3D utilities (flattened; no submodules inside the file)
    include("VEMGeo/triangulation.jl")
    include("VEMGeo/integration.jl")
    include("VEMGeo/neighbour_search.jl")

    # include("mesh.jl")
    # include("face_projector.jl")
    # include("volume_projector.jl")

    # Exports: curated public API
    export 
        # Types
        Node, NManifold, Topology, Edge, Area, Volume, RootIterator,
        Monomial, Polynomial, PolynomialBase, BaseInfo, StretchedMatrix,
        FaceTriangulations3D, VolumeIntegralData, FaceData,
        
        # Topology core
        get_coords, get_id, is_active, is_root,
        get_nodes, get_edges, get_areas, get_volumes,
        get_volume_node_ids, get_volume_area_ids,
        get_area_node_ids, get_edge_node_ids,
        add_node!, add_area!, add_volume!,
        num_roots,
        get_num_active_nodes, get_num_active_root_edges, get_num_active_root_areas, get_num_active_root_volumes,

        # Utilities
        get_next_idx, get_prev_idx, get_unique_values, find_single_intersec, max_node_distance,
        # geometry_to_vtk,

        # Lagrange utils (selective)
        lagrange_shape_function, interpolate_edge,

        # Polynomials
        derivative, ∇, ∂,
        get_base, sol_proj, grad_sol_proj, get_ϕi, ∇p, ∇x, poly_grad,

        # Stretched matrices
        stretch,

        # Triangulation / geometry helpers
        triangulate_polygon, triangulate_planar_polygon3D,
        tetrahedralize_points_convex, tetrahedralize_volume_local_ids,
        tetrahedralize_volume, build_tet_topology_from_volume,
        get_area_triangles, polygon_area3D,
        integrate_polynomial_over_face, integrate_polynomial_over_volume,
        volume_from_faces, volume_of_topo_volume,

        # Polynomial transforms
        compute_transformation_coeffs2d_to_3d, compute_transformation_coeffs3d_to_2d,
        compute_transformation_coeffs3d_to_2d!, compute_transformation_coeffs2d_to_2d!,
        compute_transformation_coeffs2d_to_2d,

        # Integration precomputes
        precompute_face_monomials, compute_face_integral, compute_face_integral_unshifted,
        get_area, get_bc, get_hf, get_outward_normal,
        precompute_volume_monomials, compute_volume_integral_unshifted,
        
        # Topology transforms
        transform_topology_planar!#, create_element_neighbour_list
end # module VEMGeo

# Re-export the VEMGeo API from the top-level module
@reexport using .VEMGeo

module VEMUtils
    using StaticArrays
    using Statistics
    using LinearAlgebra
    using FixedSizeArrays
    using ..VEMGeo
    import ..VEMGeo: FlattenVecs, iterate_volume_areas, iterate_element_edges, lagrange_shape_function, interpolate_edge
    import ..VEMGeo: get_iterative_area_vertex_ids, getexp, unit_sol_proj, stretch
    import ..VEMGeo: get_gauss_legendre, get_gauss_lobatto, base_eval
    import ..MaterialLaws: eval_hessian, eval_gradient_and_hessian, Helmholtz
    using Tensorial: ⊡₂
    using Octavian
    using Bumper
    using SparseArrays
    using WriteVTK
    import Static

    include("mesh.jl")
    include("face_projector.jl")
    include("element_mapping.jl")
    include("volume_projector.jl")
    include("VEMUtils/dof_handler.jl")
    include("VEMUtils/cell_values.jl")
    include("VEMUtils/constraint_handler.jl")
    include("VEMUtils/apply.jl")
    include("VEMUtils/assembler.jl")
    include("VEMUtils/mesh_utils.jl")
    include("VEMUtils/vtkexports.jl")
    include("cached_arrays.jl")
    include("VEMUtils/element_routine.jl")
    include("VEMUtils/post_processing.jl")
    export 
        # Mesh and element types
        Mesh, ElType, StandardEl, SerendipityEl,
        # Mesh accessors
        get_vertices, get_gauss_nodes, get_face_moments, get_volume_moments,
        add_node_set!, add_edge_set!, add_face_set!,
        # Face/volume projectors
        create_B_mat, create_D_mat, h1_projectors!, 

        # Meshing utilities
        create_rectangular_mesh, create_unit_rectangular_mesh,
        create_voronoi_mesh, extrude_to_3d, load_vtk_mesh,
        create_voronoi_mesh_3d, write_voronoi3d_vtk,
        relax_voronoi3d_seeds, relax_voronoi3d,
        # Local mapping
        NodeID2LocalID,
        # Volume assembly helpers
        create_volume_bmat, create_volume_vem_projectors,
        # Dof handler
        DofHandler, get_dof, get_dofs!, get_dofs,
        # Cell values
        CellValues, reinit!,get_n_cell_dofs,ElementBaseFunctions,
        # Constraint handler
        ConstraintHandler, add_dirichlet_bc!, add_neumann_bc!,
        # Apply
        apply!,
        # Assembler
        Assembler, local_assembly!, assemble!,
        # VTK exports
        write_vtk,
        # Cached arrays
        CachedArray, CachedMatrix, CachedVector, DefCachedMatrix, DefCachedVector, setsize!,

        # Element routine
        build_kel_quad_point!, symmetrize_kel!, build_local_kel_and_f!,
        # Post processing
        compute_error
end # module VEMUtils

@reexport using .VEMUtils

using PrecompileTools: @compile_workload    # this is a small dependency

@compile_workload begin

    include("../tests/full_problems/poisson_prec.jl")
    include("../tests/full_problems/poisson_voronoi_second_order_prec.jl")

end



end # module Ju3VEM
