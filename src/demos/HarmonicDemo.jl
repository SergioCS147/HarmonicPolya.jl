using LinearAlgebra, DynamicPolynomials, Plots
include("../HarmonicBases.jl")
using .HarmonicBases

@polyvar x[1:3]

norm = sum(x.*x)

p_1 = x[1]^2
decomp_1 = HarmonicBases.harmonicdecomposition(p_1,x)
hp_1 = sum(map(v->norm^(floor(Int,v[1]/2)) * v[2],decomp_1))

p_20 = x[1]^20
decomp_20 = HarmonicBases.harmonicdecomposition(p_20,x)
hp_20 = sum(map(v->norm^(floor(Int,v[1]/2)) * v[2],decomp_20))