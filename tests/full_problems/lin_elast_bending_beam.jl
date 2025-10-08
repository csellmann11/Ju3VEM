using StaticArrays
using Test
using Ju3VEM
using Bumper
using Random
using BenchmarkTools
using LinearAlgebra
using Ju3VEM.VEMGeo: transform_topology_planar!,transform_topology_linear_elements!
using Ju3VEM.VEMGeo: _refine!,_coarsen!,stretch
using Ju3VEM.VEMUtils: write_vtk
using Ju3VEM.VEMUtils: add_face_set!
using Ju3VEM.MaterialLaws: E_ν_to_lame, lame_to_E_ν, Ψlin

include("../topo_tests/ana_error_compute.jl")
# using Ju3VEM.VEMUtils: create_volume_vem_projectors, reinit!,get_n_cell_dofs
# using Ju3VEM.VEMUtils: add_node_set!,add_dirichlet_bc!,apply!
# using Ju3VEM.VEMUtils: Mesh, StandardEl, create_volume_bmat, h1_projectors!, create_node_mapping

const U = 3


# Grid parameters
nx,ny,nz =  (3,3,10) .*3
l_beam = 5.0
mesh2d = create_voronoi_mesh((0.0,0.0),(1.0,1.0),nx,ny,StandardEl{1},false)
mesh   = extrude_to_3d(nz,mesh2d,l_beam);
# mesh = create_unit_rectangular_mesh(nx,ny,nz, StandardEl{1}) 
# mesh = create_rectangular_mesh(
#     nx,ny,nz,
#     1.0,1.0,l_beam,StandardEl{1}
# )





function is_dirichlet_boundary(x)
    return x[3] ≈ 0  
end

function is_neumann_boundary(x,pars)
    l_beam = pars[1]
    return x[3] ≈ l_beam
end

# function is_dirichlet_boundary(x)
#     return x[1] == 0  || x[1] == l_beam
# end

# function is_neumann_boundary(x,pars)
#     l_beam = pars[1]
#     return x[3] == 1.0 && (2.0 ≤ x[1] ≤ 3.0)
# end

add_face_set!(mesh, "neumann", x -> is_neumann_boundary(x,l_beam))


add_node_set!(mesh, "dirichlet", is_dirichlet_boundary)
ch = ConstraintHandler{U}(mesh)

cv = CellValues{U}(mesh)
add_dirichlet_bc!(ch,cv.dh,"dirichlet",x -> zero(SVector{U,Float64}))



traction_force(x) = SA[0.0,1.0,0.0]
add_neumann_bc!(ch,cv.dh,cv.facedata_col,"neumann",traction_force)



function build_kel!(
    kelement::CachedMatrix{Float64},
    ebf::ElementBaseFunctions{D,O,U},
    hvol::Float64,
    bc::StaticVector{3,Float64},
    abs_volume::Float64,
    mat_law::Helmholtz) where {D,O,U}

    ∇u_dummy = zero(SMatrix{U,D,Float64,U*D})
    ℂ = eval_hessian(mat_law,∇u_dummy)
    setsize!(kelement,(length(ebf),length(ebf)))
    @no_escape begin
        grad_vals = @alloc(SMatrix{U,D,Float64,U*D},length(ebf))
        for (i,p_i) in enumerate(ebf)
            grad_vals[i] = ∇x(p_i,hvol,zero(bc))
        end

        @inbounds for (i,∇pix) in enumerate(grad_vals)
            ℂ∇pix = ℂ ⊡₂ ∇pix
            for (j,∇pjx) in enumerate(grad_vals)
                kelement.array[i,j] = abs_volume*dot(ℂ∇pix,∇pjx)
            end
        end
    end
    kelement
end

function assembly(cv::CellValues{D,U},f::F) where {D,U,F<:Function}
    ass = Assembler{Float64}(cv)

    base3d = get_base(BaseInfo{3,1,U}())
    rhs_element = DefCachedVector{Float64}()
    kelement = DefCachedMatrix{Float64}()

    E = 1e03; ν = 0.3
    λ,μ = E_ν_to_lame(E,ν)
    mat_law = Helmholtz(Ψlin,D,U,λ,μ)
    ∇u_dummy = one(SMatrix{U,D,Float64,U*D})
    γ = ∇u_dummy ⊡₂ eval_hessian(mat_law,∇u_dummy)⊡₂ ∇u_dummy


    for element in RootIterator{4}(cv.mesh.topo)
        reinit!(element.id,cv)

        hvol = cv.volume_data.hvol
        bc = cv.volume_data.vol_bc
        abs_volume = cv.volume_data.integrals[1]

        
        proj_s, proj = create_volume_vem_projectors(
            element.id,cv.mesh,cv.volume_data,cv.facedata_col,cv.vnm)


        ebf = ElementBaseFunctions(base3d,stretch(proj_s,Val(U)))
    
        kelement = build_kel!(kelement,ebf,hvol,bc,abs_volume,mat_law)

 
  
        setsize!(rhs_element,length(ebf))
        stab =  γ*(I-proj)'*(I-proj)*hvol/64
        kelement .+= stretch(stab,Val(U))



        local_assembly!(ass,kelement,rhs_element)


    end

    kglobal, rhsglobal = assemble!(ass)


  

    kglobal, rhsglobal
end


@time "assembly time" k_global,rhs_global = assembly(cv,rhs_fun);



apply!(k_global,rhs_global,ch)

println("rhs_max is $(maximum(rhs_global))")


@time "solver" u = lu(k_global) \ rhs_global;


max_idx = argmax(u)
println("max u value is $(u[max_idx])")

write_vtk(mesh.topo,"vtk/lin_elast_bending_beam", cv.dh, u)
num_elemens = length(RootIterator{4}(mesh.topo)) 
println("number of elements is $num_elemens")
println("number of nodes is $(length(get_nodes(mesh.topo)))")



