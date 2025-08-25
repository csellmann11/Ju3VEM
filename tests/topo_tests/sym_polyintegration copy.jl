using Ju3VEM
using Ju3VEM.VEMGeo: get_base_len,get_plane_parameters,face2d_sym_integral
using StaticArrays
using LinearAlgebra
using Symbolics, Chairmarks
using Test
using FixedSizeArrays, Bumper
using Ju3VEM.VEMGeo: D2FaceParametrization, project_to_2d_abs, project_to_2d_rel, project_to_3d, project_to_3d_flat


quad_nodes = [SA[0.0,0.0,0.0], SA[1.0,0.0,0.0], SA[1.0,1.0,1.0], SA[0.0,1.0,1.0]] .* 2 .+ Scalar(SA[1.0,1.0,1.0])
hf = max_node_distance(quad_nodes)
# u, v, n, nus,p0 = get_plane_parameters(quad_nodes)
plane = D2FaceParametrization(quad_nodes)



a   = face2d_sym_integral(quad_nodes,0,0,SA[0.0,0.0],plane)
bcu = face2d_sym_integral(quad_nodes,1,0,SA[0.0,0.0],plane)/a
bcv = face2d_sym_integral(quad_nodes,0,1,SA[0.0,0.0],plane)/a

face_data = precompute_face_monomials(quad_nodes,Val(5))



@test face_data.bc ≈ SA[bcu,bcv] 


x1 = quad_nodes[1]; x2 = quad_nodes[2]; x3 = quad_nodes[3]; x4 = quad_nodes[4]
x2d1 = project_to_2d_abs(x1,plane)
x2d2 = project_to_2d_abs(x2,plane)
x2d3 = project_to_2d_abs(x3,plane)
x2d4 = project_to_2d_abs(x4,plane)

quad_nodes_2d = [x2d1, x2d2, x2d3, x2d4] 

plane2d = D2FaceParametrization(quad_nodes_2d)

a2d   = face2d_sym_integral(quad_nodes_2d,0,0,SA[0.0,0.0],plane2d)
bcu2d = face2d_sym_integral(quad_nodes_2d,1,0,SA[0.0,0.0],plane2d)/a2d
bcv2d = face2d_sym_integral(quad_nodes_2d,0,1,SA[0.0,0.0],plane2d)/a2d
bc    = SA[bcu2d,bcv2d] 
bc3d = SA[2.0,2.0,2.0]


@test bcu ≈ bcu2d
@test bcv ≈ bcv2d  

bc3d = project_to_3d(face_data.bc,plane)
@test bc3d ≈ SA[2.0,2.0,2.0]

@test x2 ≈ project_to_3d(x2d2,plane)
@test x3 ≈ project_to_3d(x2d3,plane)
@test x4 ≈ project_to_3d(x2d4,plane)

true_area = (x2d2[1] - x2d1[1])*(x2d3[2] - x2d1[2])


m3d = Monomial(1.0,SA[2,1,1])


bc3d = project_to_3d(face_data.bc,face_data.plane)

bc3d_plane = project_to_3d_flat(face_data.bc,face_data.plane)


I = compute_face_integral(m3d,face_data)

b = @b compute_face_integral($m3d,$face_data)
display(b)

ref = 4*sqrt(2)/(9*hf^4)

@test I ≈ ref






quad_nodes2 = [SA[0.0,0.0,0.0], SA[1.0,0.0,0.0], SA[1.0,1.0,0.0], SA[0.0,1.0,0.0]] .* 2 .+ Scalar(SA[1.0,1.0,1.0])
hf = max_node_distance(quad_nodes2)
u, v, n, nus,p0 = Ju3VEM.VEMGeo.get_plane_parameters(quad_nodes2)

# p0 = 0*p0
x1d2 = SA[dot(quad_nodes2[1],u),dot(quad_nodes2[1],v)]
x2d2 = SA[dot(quad_nodes2[2],u),dot(quad_nodes2[2],v)]
x2d3 = SA[dot(quad_nodes2[3],u),dot(quad_nodes2[3],v)]
x2d4 = SA[dot(quad_nodes2[4],u),dot(quad_nodes2[4],v)]


_x1d2 = SA[dot(quad_nodes2[1]-p0,u),dot(quad_nodes2[1]-p0,v)]
_x2d2 = SA[dot(quad_nodes2[2]-p0,u),dot(quad_nodes2[2]-p0,v)]
_x2d3 = SA[dot(quad_nodes2[3]-p0,u),dot(quad_nodes2[3]-p0,v)]
_x2d4 = SA[dot(quad_nodes2[4]-p0,u),dot(quad_nodes2[4]-p0,v)]

p02d = SA[dot(p0,u),dot(p0,v)]

@test x1d2 == _x1d2 + p02d
@test x2d2 == _x2d2 + p02d
@test x2d3 == _x2d3 + p02d
@test x2d4 == _x2d4 + p02d




face_data2 = precompute_face_monomials(quad_nodes2,Val(5))
bc3d = project_to_3d(face_data2.bc,face_data2.plane)
bc3d_plane  = project_to_3d_flat(face_data2.bc,face_data2.plane)


m3d = Monomial(1.0,SA[0,1,1])

I = compute_face_integral(m3d,face_data2,zero(bc3d),2.0)

b = @b compute_face_integral($m3d,$face_data2,zero($bc3d),2.0)
display(b)



@test I ≈ 2.0

quad_nodes3 = [SA[0.0,0.0,0.0], SA[1.0,0.0,0.0], SA[1.0,1.0,0.0], SA[0.0,1.0,0.0]] .* 2 .+ 2*Scalar(SA[1.0,1.0,0.0])
face_data3 = precompute_face_monomials(quad_nodes3,Val(1))
b = @b precompute_face_monomials($quad_nodes3,Val(1))
display(b)



m3d = Monomial(1.0,SA[0,1,0])
m2d = Monomial(1.0,SA[0,1])

I  = compute_face_integral(m3d,face_data3,zero(bc3d),2.0)
I2 = compute_face_integral(m2d,face_data3,zero(p02d),2.0)

b = @b compute_face_integral($m2d,$face_data3,zero($p02d),2.0)
display(b)

@test I  ≈ 6.0 
@test I2 ≈ 6.0 


I3 = compute_face_integral(m2d,face_data3,face_data3.bc,2.0)
I4 = compute_face_integral_unshifted(m2d,face_data3)

@b compute_face_integral_unshifted($m2d,$face_data3)






## Testing plane projector 
using Ju3VEM.VEMGeo: D2FaceParametrization, project_to_2d_abs, project_to_2d_rel, project_to_3d, project_to_3d_flat

plane = D2FaceParametrization(quad_nodes)

x2d = project_to_2d_abs(quad_nodes[1],plane)
x3d = project_to_3d(x2d,plane)
x3d_flat = project_to_3d_flat(x2d,plane)

plane2 = D2FaceParametrization(quad_nodes2)

x2d2 = project_to_2d_abs(quad_nodes2[1],plane2)
x3d2 = project_to_3d(x2d2,plane2)
x3d_flat2 = project_to_3d_flat(x2d2,plane2)


