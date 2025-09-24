using StaticArrays
using Test
using Ju3VEM
using LinearAlgebra
using Ju3VEM.VEMGeo: stretch, unit_sol_proj
using Chairmarks
using FixedSizeArrays
using Bumper
using Ju3VEM.VEMUtils: DefCachedMatrix,DefCachedVector, setsize!, CachedMatrix, CachedVector
# using Ju3VEM.VEMUtils: add_node_set!,add_dirichlet_bc!,apply!
# using Ju3VEM.VEMUtils: Mesh, StandardEl, create_volume_bmat, h1_projectors!, create_node_mapping
function rhs_fun(x) 
    return 3*π^2*sin(π*x[1])*sin(π*x[2])*sin(π*x[3])
end

# Grid parameters
# nx = 30; ny = 30; nz = 30
nx,ny,nz = 10,10,10
mesh = create_unit_rectangular_mesh(nx,ny,nz, StandardEl{1})

add_node_set!(mesh, "dirichlet", x -> x[1] == 0 || x[1] == 1 ||
    x[2] == 0 || x[2] == 1 || x[3] == 0 || x[3] == 1)
ch = ConstraintHandler{1}(mesh)

@time "CellValues time" cv = CellValues{1}(mesh);
add_dirichlet_bc!(ch,cv.dh,"dirichlet",x -> 0.0)


function inner_prod(v1,v2,cv::CellValues)

    @assert length(v1) == length(v2) "vectors must have the same length"
    vd = cv.volume_data
    hvol = vd.hvol

    return sum(
        compute_volume_integral_unshifted(v1[i]*v2[i],vd,hvol) 
        for i in eachindex(v1,v2)
    )
end




function build_kel!(
    kelement::CachedMatrix{Float64},
    ebf::ElementBaseFunctions,
    hvol::Float64,
    bc::StaticVector{3,Float64},
    abs_volume::Float64
)

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



# @time "assembly time" k_global,rhs_global = let cv = cv, mesh = mesh
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
        for (node_id,i) in cv.vnm.map
            rhs_element[i] = f(cv.mesh[node_id])/length(ebf)*abs_volume 
        end

        stab = (I-proj)'*(I-proj)*hvol/8
        kelement .+= stab

        local_assembly!(ass,kelement,rhs_element)


    end
    kglobal, rhsglobal = assemble!(ass)


  

    kglobal, rhsglobal
end

@time "assembly time" k_global,rhs_global = assembly(cv,rhs_fun)


apply!(k_global,rhs_global,ch)
@time "solver_time" u = k_global \ rhs_global


max_idx = argmax(u)
if iseven(nx) && iseven(ny) && iseven(nz)
    @test mesh[max_idx] ≈ [0.5,0.5,0.5]
    @test u[max_idx] ≈ 1.0 atol = 0.1
end


# @time write_vtk(mesh.topo, "vtk/simple_poisson_test",cv.dh,u)