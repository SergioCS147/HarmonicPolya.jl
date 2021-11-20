module FawziFangOperator;

include("HarmonicBases.jl")
include("SphericalQuadrature.jl")
using LinearAlgebra, DynamicPolynomials, FastGaussQuadrature, PolynomialBases, SpecialFunctions, .HarmonicBases, .SphericalQuadrature

function normalizationfactor(n,deg)
    if n < 1
        throw(DomainError(n,"dimension must be positive"));
    end
    if deg < 0
        throw(DomainError(deg,"degree must be nonnegative"));
    end
    factor = pi*(2.0^(3-n))*factorial(n-3)/(SpecialFunctions.gamma((n-2)/2)^2 * (n-2)/2);
    if deg > 0
        factor = normalizationfactor(n,deg-1)*((deg-1+n-2)*(deg-1+(n-2)/2))/(deg*(deg+(n-2)/2));
    end
    return factor;
end

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
    z,wz = FastGaussQuadrature.gaussjacobi(4*m+deg,(n-3)/2,(n-3)/2);
    ωn = SphericalQuadrature.surfaceareasphere(n-1)
    ωn1 = SphericalQuadrature.surfaceareasphere(n-2)
    indices = collect(Base.product(0:1:m, 0:1:m));
    #ω1 = SphericalQuadrature.surfaceareasphere(n-1);
    #dimH(n1,j) = binomial(n1+j-1,n1-1) - binomial(n1+j-3,n1-1)
    #ω = SphericalQuadrature.surfaceareasphere(n);
    #normfactor(deg)=1/sqrt(pi*(2.0^(1-(n-2)))*SpecialFunctions.gamma(deg+(n-2))/(factorial(deg)*(deg+(n-2)/2)*SpecialFunctions.gamma((n-2)/2)^2))
    χ = ((x,i,j) -> ((1/(sqrt(normalizationfactor(n,i))*sqrt(normalizationfactor(n,j))))*PolynomialBases.gegenbauer(x,i,(n-2)/2)*PolynomialBases.gegenbauer(x,j,(n-2)/2)*h(x)) );
    return map(y->( dot(wz,χ.(z,y[1],y[2]))),indices);
end

function eigenvalues(n,m,k)
    hfun = (x -> (1/k)*sum([((1/(normalizationfactor(n,2*i)))*PolynomialBases.gegenbauer(x,2*i,n/2-1)) for i in 1:1:k]));
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

#function kernel(n,m,k)
#    hfun = (x -> (1/k)*sum([((1/(normalizationfactor(n,2*i)))*PolynomialBases.gegenbauer(x,2*i,n/2-1)) for i in 1:1:k]));
#    @polyvar s
#    h = hfun(s)
#    τh = toeplizmatrix(n,m,k,h);
#    V = eigen(τh, sortby=abs);
#    e = V.vectors[:,m+1]/norm(V.vectors[:,m+1]);
#    mons = [s^i for i in 0:1:m];
#    Φ = (dot(e,mons))^2;
#    return coefficients(Φ,[s^i for i in 0:1:(2*m)]);
#end

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
    return sum( eigvals[u[1]+1]*norm^(floor(Int,u[1]/2))*u[2] for u in desc )
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
    eigvals = eigenvalues(n,m,2*k)
    norm = sum(vars .* vars)
    return sum( (1/eigvals[u[1]+1])*norm^(floor(Int,u[1]/2))*u[2] for u in desc )
end

end

using .FawziFangOperator, DynamicPolynomials, PolynomialBases
@polyvar x[1:4]

one = x[1]^0
FawziFangOperator.toeplizmatrix(6,20,2,one)
