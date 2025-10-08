using StaticArrays
using Test
using Ju3VEM
using Bumper
using Random
using BenchmarkTools
using LinearAlgebra
using Ju3VEM.VEMGeo: _refine!,_coarsen!

include("../topo_tests/ana_error_compute.jl")
# using Ju3VEM.VEMUtils: create_volume_vem_projectors, reinit!,get_n_cell_dofs
# using Ju3VEM.VEMUtils: add_node_set!,add_dirichlet_bc!,apply!
# using Ju3VEM.VEMUtils: Mesh, StandardEl, create_volume_bmat, h1_projectors!, create_node_mapping

const U = 1

function rhs_fun(x) 
    val = 3*π^2*sin(π*x[1])*sin(π*x[2])*sin(π*x[3])
    return SVector{U}(val for i in 1:U)
    # return 0.0
end

function ana_u_poisson(x) 
    val = sin(π*x[1])*sin(π*x[2])*sin(π*x[3])
    return SVector{U}(val for i in 1:U)
    # return 1.0 + x[1] + x[2] + x[3]
end

# Grid parameters
nx,ny,nz =  8,8,8
# mesh2d = create_voronoi_mesh((0.0,0.0),(1.0,1.0),nx,ny,StandardEl{1},false)
# mesh   = extrude_to_3d(nz,mesh2d,1.0);
# mesh = create_unit_rectangular_mesh(nx,ny,nz, StandardEl{1}) 
mesh = load_vtk_mesh("tests/full_problems/voronoi4000_3d_with_lloyd.vtk", StandardEl{1});
# mesh = load_vtk_mesh("tests/full_problems/voronoi4000_3d_no_lloyd.vtk", StandardEl{1});

function rand_refinement(mesh,num_ref = 2)
    rng = Random.MersenneTwister(123)
    for i in 1:num_ref
        for vol in RootIterator{4}(mesh.topo)
            if rand(rng,0:1) |> Bool
                _refine!(vol,mesh.topo)
            end
        end
    end
end



mesh = Mesh(mesh.topo, StandardEl{1}())


function is_dirichlet_boundary(x)
    return x[1] ≈ 0 || x[1] ≈ 1 ||
                x[2] ≈ 0 || x[2] ≈ 1 || x[3] ≈ 0# || x[3] ≈ 1
end

function is_neumann_boundary(x)
    return x[3] ≈ 1 
end

add_face_set!(mesh, "neumann", is_neumann_boundary)


add_node_set!(mesh, "dirichlet", is_dirichlet_boundary)
ch = ConstraintHandler{U}(mesh)

cv = CellValues{U}(mesh)
add_dirichlet_bc!(ch,cv.dh,"dirichlet",x -> ana_u_poisson(x))




grad_u(x) = SA[
    pi*cos(π*x[1])*sin(π*x[2])*sin(π*x[3]),
    pi*sin(π*x[1])*cos(π*x[2])*sin(π*x[3]),
    pi*sin(π*x[1])*sin(π*x[2])*cos(π*x[3])
]


grad_times_n(x) = SVector{U}(grad_u(x)[3] for i in 1:U)
@time add_neumann_bc!(ch,cv.dh,cv.facedata_col,"neumann",grad_times_n)


function build_kel!(
    kelement::CachedMatrix{Float64},
    ebf::ElementBaseFunctions{D,O,U},
    hvol::Float64,
    bc::StaticVector{3,Float64},
    abs_volume::Float64) where {D,O,U}

    # kelement = FixedSizeMatrix{Float64}(undef,length(ebf),length(ebf))
    setsize!(kelement,(length(ebf),length(ebf)))
    @no_escape begin
        grad_vals = @alloc(SMatrix{U,D,Float64,U*D},length(ebf))
        for (i,p_i) in enumerate(ebf)
            grad_vals[i] = ∇x(p_i,hvol,zero(bc))
        end

        @inbounds for (i,∇pix) in enumerate(grad_vals)
            for (j,∇pjx) in enumerate(grad_vals)
                kelement.array[i,j] = abs_volume*dot(∇pix,∇pjx)
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


    for element in RootIterator{4}(cv.mesh.topo)
        reinit!(element.id,cv)

        hvol = cv.volume_data.hvol
        bc = cv.volume_data.vol_bc
        abs_volume = cv.volume_data.integrals[1]

        
        proj_s, proj = create_volume_vem_projectors(
            element.id,cv.mesh,cv.volume_data,cv.facedata_col,cv.vnm)


        ebf = ElementBaseFunctions(base3d,stretch(proj_s,Val(U)))
    
        kelement = build_kel!(kelement,ebf,hvol,bc,abs_volume)

 
        setsize!(rhs_element,length(ebf))
        fmean = f(bc)
        for (i,base_fun) in enumerate(ebf)
            rhs_element.array[i] = fmean ⋅ compute_volume_integral_unshifted(
                base_fun,cv.volume_data,hvol
            )
        end
        stab = (I-proj)'*(I-proj)*hvol/4
        kelement .+= stretch(stab,Val(U))



        local_assembly!(ass,kelement,rhs_element)


    end

    kglobal, rhsglobal = assemble!(ass)


  

    kglobal, rhsglobal
end


@time "assembly time" k_global,rhs_global = assembly(cv,rhs_fun);



apply!(k_global,rhs_global,ch)



@time "solver" u = lu(k_global) \ rhs_global;


max_idx = argmax(u)
println("max u value is $(u[max_idx])")

# write_vtk(mesh.topo,"vtk/simple_poisson_test_ref", cv.dh, u)
num_elemens = length(RootIterator{4}(mesh.topo)) 
println("number of elements is $num_elemens")
println("number of nodes is $(length(get_nodes(mesh.topo)))")
@time l2_error = compute_error(ana_u_poisson,cv,u)
println("l2 error is $l2_error")

check_mesh(cv)