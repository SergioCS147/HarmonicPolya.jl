module FawziFangOperator;

include("HarmonicBases.jl")
include("SphericalQuadrature.jl")
using LinearAlgebra, DynamicPolynomials, FastGaussQuadrature, PolynomialBases, SpecialFunctions, .HarmonicBases, .SphericalQuadrature
using Convex, SCS

function risingfactorial(a,b)
    res = 1
    if b > 0
        res = a*risingfactorial(a+1,b-1)
    elseif b < 0
        throw(DomainError(b,"b must be positive"))
    end
    return res
end

Acoefficient(n,i,j) = 2.0^((n-2)/2 + 1) * binomial(i+n-2-1,i) * binomial(j+n-2-1,j)

bcoefficient(n,i,j,r) = sum([ sum([ (r+(n-2)/2)*(-1.0)^(p+r) * (binomial(i,p-m)*binomial(j,m)*(n-2+j)^m * risingfactorial(n-2+n-j,p-m)*risingfactorial(p,r) * risingfactorial((n-2)/2 + 1/2, p) )/(risingfactorial((n-2)/2 + 1/2, m)*risingfactorial((n-2)/2 + 1/2, p-m)*risingfactorial(n-2,r+p+1)) for m in 0:1:p]) for p in r:1:n])

function eigenvalues(n,m,k)
    e = Convex.Variable(m+1)
end 

function evaluation(m,p,vars)
    desc = HarmonicBases.harmonicdecomposition(p,vars)
    n = length(vars)
    k = floor(Int,maxdegree(p)/2)
    if 2*k != maxdegree(p)
        throw(DomainError(p,"polynomial p must have even degree"))
    end
    if 2*k != mindegree(p)
        throw(DomainError(p,"polynomial p must be homogeneous"))
    end
    eigvals = eigenvalues(n,m,2*k)
    norm = sum(vars .* vars)
    return sum( eigvals[floor(Int,u[1]/2)+1]*norm^(floor(Int,u[1]/2))*u[2] for u in desc )
end

function inverseeval(m,p,vars)
    desc = HarmonicBases.harmonicdecomposition(p,vars)
    n = length(vars)
    k = floor(Int,maxdegree(p)/2)
    if 2*k != maxdegree(p)
        throw(DomainError(p,"polynomial p be of even degree"))
    end
    if 2*k != mindegree(p)
        throw(DomainError(p,"polynomial p must be homogeneous"))
    end
    eigvals = eigenvalues(n,m,k)
    norm = sum(vars .* vars)
    return sum( (1/eigvals[floor(Int,u[1]/2)+1])*norm^(floor(Int,u[1]/2))*u[2] for u in desc )
end

end