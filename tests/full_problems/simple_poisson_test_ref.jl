using StaticArrays
using Test
using Ju3VEM
using LinearAlgebra
using Ju3VEM.VEMGeo: transform_topology_planar!,transform_topology_linear_elements!
using Ju3VEM.VEMGeo: _refine!,_coarsen!
using Ju3VEM.VEMUtils: write_vtk
# using Ju3VEM.VEMUtils: create_volume_vem_projectors, reinit!,get_n_cell_dofs
# using Ju3VEM.VEMUtils: add_node_set!,add_dirichlet_bc!,apply!
# using Ju3VEM.VEMUtils: Mesh, StandardEl, create_volume_bmat, h1_projectors!, create_node_mapping
function rhs_fun(x) 
    return 3*π^2*sin(π*x[1])*sin(π*x[2])*sin(π*x[3])
end

# Grid parameters
nx,ny,nz =  10,10,10

mesh2d = create_voronoi_mesh((0.0,0.0),(1.0,1.0),nx,ny,StandardEl{1},false)
mesh   = extrude_to_3d(nz,mesh2d,1.0);
# mesh = create_unit_rectangular_mesh(nx,ny,nz, StandardEl{1})

n_element_old = length(RootIterator{4}(mesh.topo))

@time "refinement" for i in 1:2
    for vol in RootIterator{4}(mesh.topo)
        if rand(0:1) |> Bool
            _refine!(vol,mesh.topo)
        end
    end
end





mesh = Mesh(mesh.topo, StandardEl{1}())


function is_boundary(x)
    return x[1] == 0 || x[1] == 1 ||
    x[2] == 0 || x[2] == 1 || x[3] == 0 || x[3] == 1
end
add_node_set!(mesh, "dirichlet", is_boundary)
ch = ConstraintHandler{1}(mesh)

cv = CellValues{1}(mesh)
add_dirichlet_bc!(ch,cv.dh,"dirichlet",x -> 0.0)

function build_kel!(
    kelement::CachedMatrix{Float64},
    ebf::ElementBaseFunctions,
    hvol::Float64,
    bc::StaticVector{3,Float64},
    abs_volume::Float64)

    # kelement = FixedSizeMatrix{Float64}(undef,length(ebf),length(ebf))
    setsize!(kelement,(length(ebf),length(ebf)))
    @no_escape begin
        grad_vals = @alloc(SVector{3,Float64},length(ebf))
        for (i,p_i) in enumerate(ebf)
            grad_vals[i] = ∇x(p_i,hvol,zero(bc))
        end

        @inbounds for (i,∇pix) in enumerate(grad_vals)
            for (j,∇pjx) in enumerate(grad_vals)
                kelement[i,j] = abs_volume*dot(∇pix,∇pjx)
            end
        end
    end
    kelement
end

function assembly(cv::CellValues{D,U},f::F) where {D,U,F<:Function}
    ass = Assembler{Float64}(cv)

    base3d = get_base(BaseInfo{3,1,1}())
    rhs_element = DefCachedVector{Float64}()
    kelement = DefCachedMatrix{Float64}()

    for element in RootIterator{4}(cv.mesh.topo)
        reinit!(element.id,cv)

        hvol = cv.volume_data.hvol
        bc = cv.volume_data.vol_bc
        abs_volume = cv.volume_data.integrals[1]

        
        proj_s, proj = create_volume_vem_projectors(
            element.id,cv.mesh,cv.volume_data,cv.facedata_col,cv.vnm)

        ebf = ElementBaseFunctions(base3d,stretch(proj_s))
    
        kelement = build_kel!(kelement,ebf,hvol,bc,abs_volume)

 
        setsize!(rhs_element,length(ebf))
        fmean = f(bc)
        for (i,base_fun) in enumerate(ebf)
            rhs_element[i] = fmean*compute_volume_integral_unshifted(
                base_fun,cv.volume_data,hvol
            )
        end

        # for (node_id,i) in cv.vnm.map
        #     rhs_element[i] = f(cv.mesh[node_id])/length(ebf)*abs_volume 
        # end

        stab = (I-proj)'*(I-proj)*hvol/2
        kelement .+= stab

        local_assembly!(ass,kelement,rhs_element)


    end
    kglobal, rhsglobal = assemble!(ass)


  

    kglobal, rhsglobal
end

@time "assembly time" k_global,rhs_global = assembly(cv,rhs_fun)



apply!(k_global,rhs_global,ch)
@time "solver" u = k_global \ rhs_global



max_idx = argmax(u)
println("max u value is $(u[max_idx])")

write_vtk(mesh.topo,"vtk/simple_poisson_test_ref", cv.dh, u)