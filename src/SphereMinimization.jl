module SphereMinimization;

include("SphericalQuadrature.jl")
include("HarmonicBases.jl")
include("GammaOperator.jl")
using LinearAlgebra, FixedPolynomials, .SphericalQuadrature, .HarmonicBasis, .GammaOperator

function solvehomogeneous(p,vars,m)
    n = length(vars);
    deg = maxdegree(p);
    h2m = GammaOperator.inverseeval(m,p,vars);
    fixedh2m = FixedPolynomials.Polynomial{Float64}(h2m);
    z,wz = SphericalQuadrature.sphericalquadrature(n,deg+2*m);
    evals = map(x->FixedPolynomials.evaluate(fixedh2m,x),z);
    alpha = minimum(evals);
    return alpha
end

end