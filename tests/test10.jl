using LinearAlgebra, StaticArrays
using Ju3VEM
using Ju3VEM.VEMGeo: D2FaceParametrization
A = @SVector [0.0, 0.0, 0.0]
B = @SVector [2.0, 0.0, 0.0]
C = @SVector [1.5, 1.0, 0.0]
D = @SVector [0.0, 1.2, 0.0]

o2 = @SVector [0.0, 0.0, -0.1]  # offset applied to edge CD

Ap = A 
Bp = B 
Cp = C + o2 
Dp = D + o2 

quad_nodes = [Ap,Bp,Cp,Dp]

plane = D2FaceParametrization(quad_nodes)

quad_nodes2d = project_to_2d_abs.(quad_nodes,Ref(plane))
quad_nodes_recomp = project_to_3d.(quad_nodes2d,Ref(plane))

@show quad_nodes_recomp ≈ quad_nodes
# e1 = Bp - Ap
# e2 = Cp - Ap
# n = cross(e1, e2)
# triple = dot(n, Dp - Ap)
# dist = norm(n) == 0.0 ? Inf : abs(triple) / norm(n)

# println("A' = $Ap")
# println("B' = $Bp")
# println("C' = $Cp")
# println("D' = $Dp")
# println("scalar triple = $triple")
# println("orthogonal distance of D' to plane(A',B',C') = $dist")


