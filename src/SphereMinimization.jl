function upperbound(p,vars,m)
    n = length(vars);
    deg = DynamicPolynomials.maxdegree(p);
    fixedp = FixedPolynomials.Polynomial{Float64}(p);
    z,wz = sphericalquadrature(n,deg+2*m);
    evals = map(x->FixedPolynomials.evaluate(fixedp,x),z);
    alpha = minimum(evals);
    return alpha
end

function lowerboundsquares(p,vars,m)
    n = length(vars);
    deg = DynamicPolynomials.maxdegree(p);
    h2m = squaresinverseeval(m,p,vars);
    fixedh2m = FixedPolynomials.Polynomial{Float64}(h2m);
    z,wz = sphericalquadrature(n,deg+2*m);
    evals = map(x->FixedPolynomials.evaluate(fixedh2m,x),z);
    alpha = minimum(evals);
    return alpha
end

function lowerboundfawzi(p,vars,m)
    n = length(vars);
    deg = DynamicPolynomials.maxdegree(p);
    h2m = ffinverseeval(m,p,vars);
    fixedh2m = FixedPolynomials.Polynomial{Float64}(h2m);
    z,wz = sphericalquadrature(n,deg+2*m+1);
    evals = map(x->FixedPolynomials.evaluate(fixedh2m,x),z);
    alpha = minimum(evals);
    return alpha
end