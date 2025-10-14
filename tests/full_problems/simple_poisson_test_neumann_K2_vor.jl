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
using Tensorial: contract
include("../topo_tests/ana_error_compute.jl")

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
mesh = create_unit_rectangular_mesh(1,1,1, StandardEl{K})



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


function build_kel_quad_point!(
    kelement::CachedMatrix{Float64},
    rhs_element::CachedVector{Float64},
    cv::CellValues{D,U,ET},
    fe_cv::FR.CellValues,
    Nv::AbstractVector,
    ∇Nv::AbstractVector,
    f::F,
    qpoint::Int,
    x::StaticVector{3,Float64},
    mat_law::Helmholtz,
    pars...) where {D,U,ET,F<:Function}


    hvol = cv.volume_data.hvol
    bc = cv.volume_data.vol_bc
    x_scaled = (x .- bc)/hvol

    dΩ = FR.getdetJdV(fe_cv,qpoint) 

    #INFO: for linea materials the input grad is just a dummy
    ℂ = eval_hessian(mat_law,∇Nv[1],pars...)

    @inbounds for (i,∇Ni) in enumerate(∇Nv)
        Nix  = Nv[i](x_scaled)
        ∇Nixℂ = ℂ ⊡₂ ∇Ni(x_scaled) 
        rhs_element[i] += f(x)⋅Nix*dΩ
        for j in i:length(Nv)
            ∇Njx = ∇Nv[j](x_scaled)
            kelement[i,j] += ∇Nixℂ ⋅ ∇Njx *dΩ
        end
    end
    kelement
end

function symmetrize_kel!(
    kelement::CachedMatrix{Float64})

    @assert size(kelement,1) == size(kelement,2) "kelement must be square"
    @inbounds for i in axes(kelement,1)
        for j in (i+1):size(kelement,2)  
            kelement[j,i] = kelement[i,j]  
        end
    end
    kelement
end

function build_local_kel_and_f!(
    kelement::CachedMatrix{Float64},
    rhs_element::CachedVector{Float64},
    cv::CellValues{D,U,ET},
    element_id::Int,
    f::F,
    fe_cv::FR.CellValues,
    mat_law::Helmholtz) where {D,U,F<:Function,K,ET<:ElType{K}}


    hvol = cv.volume_data.hvol

    base3d = get_base(BaseInfo{3,K,U}())
    L = length(base3d)
    poly_base3d = get_base(BaseInfo{3,K-1,U}())
    L_grad = length(poly_base3d)

    
    proj_s, proj = create_volume_vem_projectors(
        element_id,cv.mesh,cv.volume_data,cv.facedata_col,cv.vnm)

    ebf = ElementBaseFunctions(base3d,stretch(proj_s,Val(U)))

    poly_type = eltype(ebf)
    grad_type = Core.Compiler.return_type(∇p, Tuple{poly_type,Float64})

  
    ∇Nv = FixedSizeVector{grad_type}(undef,length(ebf))
    Nv  = FixedSizeVector{poly_type}(undef,length(ebf))

    for (i,p) in enumerate(ebf)
        ∇Nv[i] = ∇p(p,hvol)
        Nv[i]  = p |> poly_type
    end

    tets_local, l2g = tetrahedralize_volume(cv.mesh.topo, element_id)

    setsize!(kelement,(length(ebf),length(ebf)))
    setsize!(rhs_element,(length(ebf),))

    for tet_local_ids in tets_local
        tets_global_ids = SVector{4,Int}(l2g[id] for id in tet_local_ids)
        tet_nodes = FR.Vec{D}.(get_nodes(cv.mesh.topo)[tets_global_ids])
        FR.reinit!(fe_cv,tet_nodes)

        for qpoint in 1:FR.getnquadpoints(fe_cv)
            
            x = FR.spatial_coordinate(fe_cv,qpoint,tet_nodes) |> SVector{3,Float64}
            build_kel_quad_point!(kelement,rhs_element,cv,fe_cv,Nv,∇Nv,f,qpoint,x,mat_law)
        end 
    end

    symmetrize_kel!(kelement)

    stab       = (I-proj)'*(I-proj)*hvol/4
    kelement .+= stretch(stab,Val(U))
end


function assembly(cv::CellValues{D,U,ET},
    f::F,
    fe_cv::FR.CellValues) where {D,U,F<:Function,K,ET<:ElType{K}}
    ass = Assembler{Float64}(cv)


    mat_law = Helmholtz(Ψpoisson,D,U)
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





