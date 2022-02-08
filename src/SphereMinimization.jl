module SphereMinimization;

include("SphericalQuadrature.jl")
include("HarmonicBases.jl")
include("GammaOperator.jl")
include("FawziFangOperator.jl")
using LinearAlgebra, DynamicPolynomials, FixedPolynomials, .SphericalQuadrature, .HarmonicBases, .GammaOperator, .FawziFangOperator

function solvehomogeneousgamma(p,vars,m)
    n = length(vars);
    deg = DynamicPolynomials.maxdegree(p);
    h2m = GammaOperator.inverseeval(m,p,vars);
    fixedh2m = FixedPolynomials.Polynomial{Float64}(h2m);
    z,wz = SphericalQuadrature.sphericalquadrature(n,deg+2*m);
    evals = map(x->FixedPolynomials.evaluate(fixedh2m,x),z);
    alpha = minimum(evals);
    return alpha
end

function solvehomogeneousfawzi(p,vars,m)
    n = length(vars);
    deg = DynamicPolynomials.maxdegree(p);
    h2m = FawziFangOperator.inverseeval(m,p,vars);
    fixedh2m = FixedPolynomials.Polynomial{Float64}(h2m);
    z,wz = SphericalQuadrature.sphericalquadrature(n,deg+2*m+1);
    evals = map(x->FixedPolynomials.evaluate(fixedh2m,x),z);
    alpha = minimum(evals);
    return alpha
end

function solvehomogeneousfauzialt(p,vars,m)
    n = length(vars);
    deg = DynamicPolynomials.maxdegree(p);
    kernel = FauziFangOperator.kernel(n,m,floor(Int,deg/2));
    h2m = sum([kernel[2*i+1]*GammaOperator.inverseeval(2*i,p,vars) for i in floor(Int,deg/2):1:floor(Int,m/2)]);
    fixedh2m = FixedPolynomials.Polynomial{Float64}(h2m);
    z,wz = SphericalQuadrature.sphericalquadrature(n,deg+2*m);
    evals = map(x->FixedPolynomials.evaluate(fixedh2m,x),z);
    alpha = minimum(evals);
    return alpha
end

end