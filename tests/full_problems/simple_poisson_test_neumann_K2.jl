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
    return -2.0
    return 0.0
end

function ana_u_poisson(x) 
    val = sin(π*x[1])*sin(π*x[2])*sin(π*x[3])
    # return SVector{U}(val for i in 1:U)
    return x[3]^2
    return 1.0 + x[1] + x[2] + x[3]
end

# Grid parameters
nx,ny,nz =  (1,1,1) .* 16
mesh2d     = create_voronoi_mesh((0.0,0.0),(1.0,1.0),nx,ny,StandardEl{K},false)
mesh_old   = extrude_to_3d(nz,mesh2d,1.0);
# mesh = create_unit_rectangular_mesh(nx,ny,nz, StandardEl{K}) 


add_face_set!(mesh_old, "dirichlet", is_dirichlet_boundary)

# transform_fun(x) = SA[x[1],x[2],x[3]^2*(x[1]+x[2]) +x[3]]
function transform_fun(x)

    x3 = if x[3] <= 0.5 
        -0.5*(x[1]+x[2])+x[3]
    else
        0.5*(x[1]+x[2])+x[3]
    end
    return SA[x[1],x[2],x3]
end
for i in eachindex(get_nodes(mesh_old.topo))
    coord = transform_fun(get_nodes(mesh_old.topo)[i])
    old_node = get_nodes(mesh_old.topo)[i]
    get_nodes(mesh_old.topo)[i] = Node(old_node.id, coord)
end

mesh = Mesh(mesh_old.topo, StandardEl{K}())
mesh.face_sets["dirichlet"] = mesh_old.face_sets["dirichlet"]


#! Warning: meshes are incorrect. Some inner nodes touch only 1 element
# mesh = load_vtk_mesh("tests/full_problems/voronoi4000_3d_with_lloyd.vtk", StandardEl{K});
# mesh = load_vtk_mesh("tests/full_problems/voronoi64_3d_no_lloyd.vtk", StandardEl{K});
# mesh = load_vtk_mesh("tests/full_problems/test_512_mesh.vtk", StandardEl{K});

function is_dirichlet_boundary(x)
    return x[1] ≈ 0 || x[1] ≈ 1 ||
                x[2] ≈ 0 || x[2] ≈ 1 || x[3] ≈ 0 || x[3] ≈ 1
end


ch = ConstraintHandler{U}(mesh)

cv = CellValues{U}(mesh)
add_dirichlet_bc!(ch,cv.dh,cv.facedata_col,"dirichlet",x -> ana_u_poisson(x))


function build_kel_quad_point!(
    kelement::CachedMatrix{Float64},
    rhs_element::CachedVector{Float64},
    ebf::ElementBaseFunctions{D,O,U},
    source_fun::F,
    x::StaticVector{3,Float64},
    x_scaled::StaticVector{3,Float64},
    dΩ::Float64,
    hvol::Float64) where {D,O,U,F<:Function}



    ∇p_vals = FixedSizeVector{SMatrix{U,D,Float64,U*D}}(undef,length(ebf))
    p_vals  = FixedSizeVector{SVector{U,Float64}}(undef,length(ebf))
    for (i,p) in enumerate(ebf)
        ∇p_vals[i] = ∇x(p,hvol,x_scaled)
        p_vals[i]  = p(x_scaled) |> SVector{U,Float64}
    end

 

    @inbounds for (i,∇pix) in enumerate(∇p_vals)
        rhs_element[i] += source_fun(x)⋅p_vals[i]*dΩ
        for j in i:length(ebf)
            ∇pjx = ∇p_vals[j]
            kelement[i,j] += ∇pix ⋅ ∇pjx *dΩ
        end
    end
    kelement
end

function symmetric_build_kel!(
    kelement::CachedMatrix{Float64})
    for i in axes(kelement,1)
        for j in (i+1):size(kelement,2)  
            kelement[j,i] = kelement[i,j]  
        end
    end
    kelement
end



function assembly(cv::CellValues{D,U,ET},
    f::F,
    fe_cv::FR.CellValues) where {D,U,F<:Function,K,ET<:ElType{K}}
    ass = Assembler{Float64}(cv)

    base3d = get_base(BaseInfo{3,K,U}())
    L = length(get_base(BaseInfo{3,K-1,U}()))
    rhs_element = DefCachedVector{Float64}()
    kelement = DefCachedMatrix{Float64}()


    for element in RootIterator{4}(cv.mesh.topo)
        reinit!(element.id,cv)

        hvol = cv.volume_data.hvol
        bc = cv.volume_data.vol_bc

        
        proj_s, proj = create_volume_vem_projectors(
            element.id,cv.mesh,cv.volume_data,cv.facedata_col,cv.vnm)


  
        ebf = ElementBaseFunctions(base3d,stretch(proj_s,Val(U)))
    
        tets_local, l2g = tetrahedralize_volume(cv.mesh.topo, element.id)

        setsize!(kelement,(length(ebf),length(ebf)))
        setsize!(rhs_element,(length(ebf),))


        for (i,tet_local_ids) in enumerate(tets_local) 
            tets_global_ids = SVector{4,Int}(l2g[id] for id in tet_local_ids)
            tet_nodes = FR.Vec{D}.(get_nodes(cv.mesh.topo)[tets_global_ids])
            FR.reinit!(fe_cv,tet_nodes)

            for qpoint in 1:FR.getnquadpoints(fe_cv)
                dΩ = FR.getdetJdV(fe_cv,qpoint)
                x = FR.spatial_coordinate(fe_cv,qpoint,tet_nodes) |> SVector{3,Float64}
                x_scaled = (x .- bc)/hvol
                build_kel_quad_point!(kelement,rhs_element,ebf,f,x,x_scaled,dΩ,hvol)
            end
        end

        symmetric_build_kel!(kelement)


        stab       = (I-proj)'*(I-proj)*hvol/4
        kelement .+= stretch(stab,Val(U))
  
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

if K == 2
    @test l2_error ≈ 0.0 atol = 1e-12
end






