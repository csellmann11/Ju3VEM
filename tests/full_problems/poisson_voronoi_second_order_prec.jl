using StaticArrays
using Test
using Ju3VEM
using Bumper
using Random
using BenchmarkTools
using LinearAlgebra
using FixedSizeArrays
import Ferrite as FR
using Ju3VEM.VEMUtils: write_voronoi3d_vtk, relax_voronoi3d, relax_voronoi3d_seeds
using Ju3VEM.VEMUtils: build_local_kel_and_f!, compute_error
using Tensorial: contract
using Ju3VEM.VEMGeo: _refine!,_coarsen!

const U = 1
const K = 2

# INFO: FEM Utils for subtrias 
fe_ip = FR.Lagrange{FR.RefTetrahedron, 1}()^U
fe_qr = FR.QuadratureRule{FR.RefTetrahedron}(K)
fe_cv = FR.CellValues(fe_qr,fe_ip)

function rhs_fun_order2(x) 
    val = -10.0
    return SVector{U}(val for i in 1:U)
end

function ana_u_poisson_order2(x) 
    val = 1 + x[1]^2 + x[2]^2 + 3*x[3]^2
    return SVector{U}(val for i in 1:U)
end

n = 1   
mesh2d = create_voronoi_mesh((0.0,0.0),(1.0,1.0),n,n,StandardEl{K},true)
mesh   = extrude_to_3d(n,mesh2d,1.0);


function rand_refinement(mesh::Mesh{D,ET},num_ref = 2) where {D,ET<:ElType{K}}
    rng = Random.MersenneTwister(123)
    for _ in 1:num_ref
        for vol in RootIterator{4}(mesh.topo)
            ref_bonus = vol.refinement_level * 0.2
            if rand(rng) + ref_bonus > 0.7
                _refine!(vol,mesh.topo)
            end
        end
    end
    mesh = Mesh(mesh.topo, StandardEl{K}())
    mesh
end

mesh = rand_refinement(mesh,1)







function is_dirichlet_boundary_order2(x)
    return x[1] ≈ 0 || x[1] ≈ 1 ||
                x[2] ≈ 0 || x[2] ≈ 1 || x[3] ≈ 0 
end

function is_neumann_boundary_order2(x)
    return x[3] ≈ 1
end

add_face_set!(mesh, "dirichlet", is_dirichlet_boundary_order2)
add_face_set!(mesh, "neumann", is_neumann_boundary_order2)

ch = ConstraintHandler{U}(mesh)

cv = CellValues{U}(mesh)
add_dirichlet_bc!(ch,cv.dh,cv.facedata_col,"dirichlet",x -> ana_u_poisson_order2(x))

add_neumann_bc!(ch,cv.dh,cv.facedata_col,"neumann",x -> SVector{U}(6*x[3] for i in 1:U))



function assembly_order2(cv::CellValues{D,U,ET},
    f::F,
    fe_cv::FR.CellValues) where {D,U,F<:Function,K,ET<:ElType{K}}
    ass = Assembler{Float64}(cv)


    mat_law     = Helmholtz{U,D}(Ψpoisson)
    rhs_element = DefCachedVector{Float64}()
    kelement    = DefCachedMatrix{Float64}()


    for element in RootIterator{4}(cv.mesh.topo)
        reinit!(element.id,cv)

        build_local_kel_and_f!(kelement,rhs_element,cv,element.id,f,fe_cv,mat_law)
  
        local_assembly!(ass,kelement,rhs_element)
    end

    kglobal, rhsglobal = assemble!(ass)

    kglobal, rhsglobal
end


@time "assembly_order2 time" k_global,rhs_global = assembly_order2(cv,rhs_fun_order2,fe_cv);



apply!(k_global,rhs_global,ch)








