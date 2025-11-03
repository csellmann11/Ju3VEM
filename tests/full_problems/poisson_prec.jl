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

const U = 1
const K = 1

# INFO: FEM Utils for subtrias 
fe_ip = FR.Lagrange{FR.RefTetrahedron, 1}()^U
fe_qr = FR.QuadratureRule{FR.RefTetrahedron}(K)
fe_cv = FR.CellValues(fe_qr,fe_ip)

function rhs_fun(x) 
    val = 3π^2*prod(sinpi.(x))
    return SVector{U}(val for i in 1:U)
end

function ana_u_poisson(x) 
    val = prod(sinpi.(x))
    return SVector{U}(val for i in 1:U)
end

n = 2
mesh = create_unit_rectangular_mesh(n,n,n, StandardEl{K})



function is_dirichlet_boundary(x)
    return x[1] ≈ 0 || x[1] ≈ 1 ||
                x[2] ≈ 0 || x[2] ≈ 1 || x[3] ≈ 0 || x[3] ≈ 1
end

add_face_set!(mesh, "dirichlet", is_dirichlet_boundary)

ch = ConstraintHandler{U}(mesh)

cv = CellValues{U}(mesh)
add_dirichlet_bc!(ch,cv.dh,cv.facedata_col,"dirichlet",x -> ana_u_poisson(x))



function assembly(cv::CellValues{D,U,ET},
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


@time "assembly time" k_global,rhs_global = assembly(cv,rhs_fun,fe_cv)
apply!(k_global,rhs_global,ch)


