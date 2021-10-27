module FauziFangOperator;

include("HarmonicBases.jl")
include("SphericalQuadrature.jl")
using LinearAlgebra, DynamicPolynomials, FastGaussQuadrature, PolynomialBases, .HarmonicBases, .SphericalQuadrature

function toeplizmatrix(n,m,k,h) 
    if m < 0
        throw(DomainError(m,"order must be non-negative"));
    end
    if n < 1
        throw(DomainError(n,"dimension must be positive"));
    end
    if k > m
        throw(DomainError(k,"degree must be lesser than or equal to m"));
    end 
    deg = DynamicPolynomials.maxdegree(h)
    z,wz = FastGaussQuadrature.gaussjacobi(4*k+deg,(n-3)/2,(n-3)/2);
    ωn = SphericalQuadrature.surfaceareasphere(n-1)
    ωn1 = SphericalQuadrature.surfaceareasphere(n-2)
    indices = collect(Base.product(0:1:m, 0:1:m));
    χ = ((x,i,j) -> (ωn1/ωn)*(PolynomialBases.gegenbauer(x,i,(n-2)/2)*PolynomialBases.gegenbauer(x,j,(n-2)/2)*h(x))/sqrt(PolynomialBases.gegenbauer(1.0,i,(n-2)/2)*PolynomialBases.gegenbauer(1.0,j,(n-2)/2)));
    ω1 = SphericalQuadrature.surfaceareasphere(n-1)
    ω = SphericalQuadrature.surfaceareasphere(n)
    return map(y->((ω1/ω)*dot(wz,χ.(z,y[1],y[2]))),indices);
end

function eigenvalues(n,m,k)
    hfun = (x -> (1.0/k)*sum([PolynomialBases.gegenbauer(x,2*i,n/2-1) for i in 1:1:k]));
    @polyvar s
    h = hfun(s)
    τh = toeplizmatrix(n,m,k,h);
    V = eigen(τh, sortby=abs);
    e = V.vectors[:,m+1]/norm(V.vectors[:,m+1]);
    gfun = (x->PolynomialBases.gegenbauer(x,0,n/2-1));
    g = gfun(s);
    values = [1.0];
    for i in 1:1:k
        gfun = (x->PolynomialBases.gegenbauer(x,2*i,n/2-1));
        g = gfun(s);
        lambda = dot(e,toeplizmatrix(n,m,k,g)*e);
        push!(values, lambda);
    end
    return values;
end

function evaluation(m,p,vars)
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

using .FauziFangOperator, DynamicPolynomials, PolynomialBases
@polyvar x[1:3]

one = x[1]^0
FauziFangOperator.toeplizmatrix(3,3,2,one)
