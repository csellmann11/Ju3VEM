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
const K = 1

# INFO: FEM Utils for subtrias 
fe_ip = FR.Lagrange{FR.RefTetrahedron, 1}()^U
fe_qr = FR.QuadratureRule{FR.RefTetrahedron}(K)
fe_cv = FR.CellValues(fe_qr,fe_ip)

function rhs_fun(x) 
    val = 0.0
    return SVector{U}(val for i in 1:U)
end

function ana_u_poisson(x) 
    val = 1 + x[1] + x[2] + x[3]
    return SVector{U}(val for i in 1:U)
end

n = 4   
mesh2d = create_voronoi_mesh((0.0,0.0),(1.0,1.0),n,n,StandardEl{K},true)
mesh   = extrude_to_3d(n,mesh2d,1.0);


# left  = (0.0, 0.0, 0.0)
# right = (1.0, 1.0, 1.0)
# rng   = Random.MersenneTwister(1)
# seeds = [SA[rand(rng), rand(rng), rand(rng)] for _ in 1:64]
# seeds = relax_voronoi3d_seeds(left, right, seeds; maxiters=15, move_tol=1e-7, step=0.8)
# mesh = create_voronoi_mesh_3d(left, right, seeds, StandardEl{K}; dedup_tol=1e-10)



function rand_refinement(mesh::Mesh{D,ET},num_ref = 2) where {D,ET<:ElType{K}}
    rng = Random.MersenneTwister(123)
    for _ in 1:num_ref
        for vol in RootIterator{4}(mesh.topo)
            ref_bonus = vol.refinement_level * 0.2
            if rand(rng) + ref_bonus > 0.7
                _refine!(vol,mesh.topo, non_planar = false)
                # Ju3VEM.VEMGeo.refine_voronoi!(vol, mesh.topo; n_subdivisions=2)  # Creates 8 child cells
            end
        end
    end
    mesh = Mesh(mesh.topo, StandardEl{K}())
    mesh
end

mesh = rand_refinement(mesh,0)







function is_dirichlet_boundary(x)
    return x[1] ≈ 0 || x[1] ≈ 1 ||
                x[2] ≈ 0 || x[2] ≈ 1 || x[3] ≈ 0 
end

function is_neumann_boundary(x)
    return x[3] ≈ 1
end

add_face_set!(mesh, "dirichlet", is_dirichlet_boundary)
add_face_set!(mesh, "neumann", is_neumann_boundary)

ch = ConstraintHandler{U}(mesh)

cv = CellValues{U}(mesh)
add_dirichlet_bc!(ch,cv.dh,cv.facedata_col,"dirichlet",x -> ana_u_poisson(x))

add_neumann_bc!(ch,cv.dh,cv.facedata_col,"neumann",x -> SVector{U}(1.0 for i in 1:U))



function assembly(cv::CellValues{D,U,ET},
    f::F,
    fe_cv::FR.CellValues) where {D,U,F<:Function,K,ET<:ElType{K}}
    ass = Assembler{Float64}(cv)


    mat_law     = Helmholtz(Ψpoisson,U,D)
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


@time "assembly time" k_global,rhs_global = assembly(cv,rhs_fun,fe_cv);



apply!(k_global,rhs_global,ch)


@time "solver" u = cholesky(Symmetric(k_global)) \ rhs_global; 


# write_vtk(mesh.topo,"vtk/poisson_voronoi_ref", cv.dh, u)
write_vtk(mesh.topo,"vtk/poisson_voronoi_ref")

l2_error,h1_error = compute_error(ana_u_poisson,cv,u)
println("l2 error is $l2_error")
println("h1 error is $h1_error")

@test h1_error ≈ 0.0 atol = 1e-12
@test l2_error ≈ 0.0 atol = 1e-12


@time "neigh_list" neigh_list = create_element_neighbour_list(mesh.topo);

hist_dict = Dict{Int,Int}()
hist_vect_dict = Dict{Int,Vector{Int}}()
for (key, value) in neigh_list
    hist_dict[length(value)] = get!(hist_dict, length(value), 0) + 1
    push!(get!(hist_vect_dict, length(value), Int[]), key)
end
println(sort(hist_dict, by = values))




