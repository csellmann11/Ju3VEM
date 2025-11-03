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
mesh = create_unit_rectangular_mesh(n,n,n, StandardEl{K})

function rand_refinement(mesh::Mesh{D,ET},num_ref = 2) where {D,ET<:ElType{K}}
    rng = Random.MersenneTwister(123)
    for _ in 1:num_ref
        elements_marked = [false for _ in 1:length(get_volumes(mesh.topo))]
        for vol in RootIterator{4}(mesh.topo)
            @assert elements_marked[vol.id] == false
            ref_bonus = vol.refinement_level * 0.2
            if rand(rng) + ref_bonus > 0.7
                _refine!(vol,mesh.topo)
            else 
                if vol.parent_id != 0 && rand(rng) - ref_bonus > 0.5
                    parent = get_volumes(mesh.topo)[vol.parent_id]
                    child_ids = parent.childs 
                    elements_marked[child_ids] .= true
                    _coarsen!(vol,mesh.topo)
                end
            end
        end
    end
    mesh = Mesh(mesh.topo, StandardEl{K}())
    mesh
end

mesh = rand_refinement(mesh,3)




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


u = cholesky(Symmetric(k_global)) \ rhs_global; 


write_vtk(mesh.topo,"vtk/simple_poisson_test_refinement", cv.dh, u)

l2_error,h1_error = compute_error(ana_u_poisson,cv,u)
println("l2 error is $l2_error")
println("h1 error is $h1_error")

@test h1_error ≈ 0.0 atol = 1e-12
@test l2_error ≈ 0.0 atol = 1e-12







