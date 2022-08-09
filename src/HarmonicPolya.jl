module HarmonicPolya;

using LinearAlgebra
using MultivariatePolynomials
using DynamicPolynomials
using FastGaussQuadrature
using SpecialFunctions
using PolynomialBases
using FixedPolynomials

const MP = MultivariatePolynomials
const DP = DynamicPolynomials

import LinearAlgebra: length

include("HarmonicBases.jl")
include("SphericalQuadrature.jl")
include("SquaresOperator.jl")
include("FawziFangOperator.jl")
include("SphereMinimization.jl")

export harmonicdecomposition
export generatebasissphere
export laplacian
export sphericalquadrature
export upperbound 
export lowerboundfawzi
export lowerboundsquares

end