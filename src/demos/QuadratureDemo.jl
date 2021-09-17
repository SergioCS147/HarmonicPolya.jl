using LinearAlgebra, FixedPolynomials, Plots
import DynamicPolynomials: @polyvar
include("../SphericalQuadrature.jl")
using .SphericalQuadrature

@polyvar x[1:4]

p_1 = Polynomial{Float64}(x[1]^4)
p_2 = Polynomial{Float64}(x[3]^3 + x[4])
p_3 = Polynomial{Float64}(x[2]^2)

z, wz = SphericalQuadrature.sphericalquadrature(4,5)

I_0 = sum(wz)
I_1 = dot(wz,p_1.(z))
I_2 = dot(wz,p_2.(z))
I_3 = dot(wz,p_3.(z))

z100,wz100 = SphericalQuadrature.sphericalquadrature(3,100)

p_4 = Polynomial{Float64}(x[1]^100)
p_5 = Polynomial{Float64}(x[1]^99)

I_4 = dot(wz100,p_4.(z100))
I_5 = dot(wz100,p_5.(z100))