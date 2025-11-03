using StaticArrays
using Test
using Ju3VEM
using Bumper
using Random
using BenchmarkTools
using LinearAlgebra
using Ju3VEM.VEMGeo: _refine!,_coarsen!
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
    val = 3*π^2*sin(π*x[1])*sin(π*x[2])*sin(π*x[3])
    # return SVector{U}(val for i in 1:U)
    return 0.0
end

function ana_u_poisson(x) 
    val = sin(π*x[1])*sin(π*x[2])*sin(π*x[3])
    # return SVector{U}(val for i in 1:U)
    return x[3]
    # return x[3]^2
end

# Grid parameters
# using Ju3VEM
left  = (0.0, 0.0, 0.0)
right = (1.0, 1.0, 1.0)
rng   = Random.MersenneTwister(1)
seeds = [SA[rand(rng), rand(rng), rand(rng)] for _ in 1:8*8]
seeds = relax_voronoi3d_seeds(left, right, seeds; maxiters=15, move_tol=1e-7, step=0.8)



# mesh  = write_voronoi3d_vtk(left, right, seeds, "vtk/my_voronoi_mesh", StandardEl{K}; dedup_tol=1e-10)
# mesh = relax_voronoi3d(left, right, seeds, StandardEl{K}; dedup_tol=1e-10)
# mesh = create_voronoi_mesh_3d(left, right, seeds, StandardEl{K}; dedup_tol=1e-10)
mesh = create_unit_rectangular_mesh(10,10,10, StandardEl{K})



function is_dirichlet_boundary(x)
    return x[1] ≈ 0 || x[1] ≈ 1 ||
                x[2] ≈ 0 || x[2] ≈ 1 || x[3] ≈ 0 #|| x[3] ≈ 1
end

function is_neumann_boundary(x)
    return x[3] ≈ 1
end


add_face_set!(mesh, "dirichlet", is_dirichlet_boundary)
add_face_set!(mesh, "neumann", is_neumann_boundary)




ch = ConstraintHandler{U}(mesh)

cv = CellValues{U}(mesh)
add_dirichlet_bc!(ch,cv.dh,cv.facedata_col,"dirichlet",x -> ana_u_poisson(x))
add_neumann_bc!(ch,cv.dh,cv.facedata_col,"neumann",x -> 1.0)



function assembly(cv::CellValues{D,U,ET},
    f::F,
    fe_cv::FR.CellValues) where {D,U,F<:Function,K,ET<:ElType{K}}
    ass = Assembler{Float64}(cv)


    mat_law = Helmholtz(Ψpoisson,U,D)
    rhs_element = DefCachedVector{Float64}()
    kelement = DefCachedMatrix{Float64}()


    for element in RootIterator{4}(cv.mesh.topo)
        reinit!(element.id,cv)

        build_local_kel_and_f!(kelement,rhs_element,cv,element.id,f,fe_cv,mat_law)
  
  
        local_assembly!(ass,kelement,rhs_element)
    end

    kglobal, rhsglobal = assemble!(ass)

    kglobal, rhsglobal
end


@time "assembly time" k_global,rhs_global = assembly(cv,rhs_fun,fe_cv);



apply!(k_global,rhs_global,ch)



@time "solver" u = cholesky(Symmetric(k_global)) \ rhs_global;
# @time "solver" u = lu(k_global) \ rhs_global;


max_idx = argmax(u)
println("max u value is $(u[max_idx])")

if K == 1
    write_vtk(mesh.topo,"vtk/poisson_K2", cv.dh, u)
end
num_elemens = length(RootIterator{4}(mesh.topo)) 
println("number of elements is $num_elemens")
println("number of nodes is $(length(get_nodes(mesh.topo)))")




@time l2_error = compute_error(ana_u_poisson,cv,u)
println("l2 error is $l2_error")





