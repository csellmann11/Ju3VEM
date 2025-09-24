using StaticArrays
using Test
using Ju3VEM
using LinearAlgebra
using Ju3VEM.VEMUtils: extrude_to_3d
function rhs_fun(x) 
    return 3*π^2*sin(π*x[1])*sin(π*x[2])*sin(π*x[3])
end

# Grid parameters
# nx = 30; ny = 30; nz = 30
nx,ny,nz =  10,10,10
dx = 1/nx; dy = 1/ny; dz = 1/nz

function is_boundary(x)
    return x[1] == 0 || x[1] == 1 ||
    x[2] == 0 || x[2] == 1 || x[3] == 0 || x[3] == 1
end


mesh2d = create_voronoi_mesh((0.0,0.0),(1.0,1.0),nx,ny,StandardEl{1},false)
mesh   = extrude_to_3d(nz,mesh2d,1.0);

add_node_set!(mesh, "dirichlet", is_boundary)
ch = ConstraintHandler{1}(mesh)

cv = CellValues{1}(mesh)
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


k_global,rhs_global = let 
    ass = Assembler{Float64}(cv)
    mesh = cv.mesh; topo = mesh.topo
    base3d = get_base(BaseInfo{3,1,1}()).base
    for element in RootIterator{4}(topo)
        reinit!(element.id,cv)


        hvol = cv.volume_data.hvol
        abs_volume = cv.volume_data.integrals[1]

        proj_s, proj = create_volume_vem_projectors(
            element.id,mesh,cv.volume_data,cv.facedata_col,cv.vnm)

        n_cell_dofs = get_n_cell_dofs(cv)

       

        
        gelement = zeros(length(base3d),length(base3d))
        for (i,m3d_i) in enumerate(base3d)
            ∇m3d_i = ∇(m3d_i,hvol)
            for (j,m3d_j) in enumerate(base3d)
                ∇m3d_j = ∇(m3d_j,hvol)
                gelement[i,j] = inner_prod(∇m3d_i,∇m3d_j,cv)
            end
        end


        kelement = proj_s' * gelement * proj_s
        rhs_element = zeros(n_cell_dofs)
        for (node_id,i) in cv.vnm.map
            x = cv.mesh[node_id]
            rhs_element[i] = rhs_fun(x)/n_cell_dofs*abs_volume 
        end

        stab = (I-proj)'*(I-proj)*hvol/8



        kelement .+= stab



        local_assembly!(ass,kelement,rhs_element)
    


    end
    kglobal, rhsglobal = assemble!(ass)


  

    kglobal, rhsglobal
end


apply!(k_global,rhs_global,ch)
@time u = k_global \ rhs_global


max_idx = argmax(u)
if iseven(nx) && iseven(ny) && iseven(nz)
    # @test mesh[max_idx] ≈ [0.5,0.5,0.5]
    # @test u[max_idx] ≈ 1.0 atol = 0.1
end


write_vtk(mesh.topo, "vtk/simple_poisson_test_distort",cv.dh,u)