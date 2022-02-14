module SphereMinimization;

include("SphericalQuadrature.jl")
include("HarmonicBases.jl")
include("SquaresOperator.jl")
include("FawziFangOperator.jl")
using LinearAlgebra, DynamicPolynomials, FixedPolynomials, .SphericalQuadrature, .HarmonicBases, .SquaresOperator, .FawziFangOperator

function lowerboundsquares(p,vars,m)
    n = length(vars);
    deg = DynamicPolynomials.maxdegree(p);
    h2m = SquaresOperator.inverseeval(m,p,vars);
    fixedh2m = FixedPolynomials.Polynomial{Float64}(h2m);
    z,wz = SphericalQuadrature.sphericalquadrature(n,deg+2*m);
    evals = map(x->FixedPolynomials.evaluate(fixedh2m,x),z);
    alpha = minimum(evals);
    return alpha
end

function lowerboundfawzi(p,vars,m)
    n = length(vars);
    deg = DynamicPolynomials.maxdegree(p);
    h2m = FawziFangOperator.inverseeval(m,p,vars);
    fixedh2m = FixedPolynomials.Polynomial{Float64}(h2m);
    z,wz = SphericalQuadrature.sphericalquadrature(n,deg+2*m+1);
    evals = map(x->FixedPolynomials.evaluate(fixedh2m,x),z);
    alpha = minimum(evals);
    return alpha
end

end