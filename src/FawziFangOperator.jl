function normalizationfactor(n,deg)
    if n < 1
        throw(DomainError(n,"dimension must be positive"));
    end
    if deg < 0
        throw(DomainError(deg,"degree must be nonnegative"));
    end
    factor = pi*(2.0^(3-n))*SpecialFunctions.gamma(n-2)/(SpecialFunctions.gamma((n-2)/2)^2 * (n-2)/2);
    if deg > 0
        factor = normalizationfactor(n,deg-1)*((deg-1+n-2)*(deg-1+(n-2)/2))/(deg*(deg+(n-2)/2));
    end
    return factor;
end

function reproducinggegenbauer(x,n,deg)
    factor = (binomial(n+deg-1,deg)-binomial(n+deg-3,deg-2))/(surfaceareasphere(n)*PolynomialBases.gegenbauer(1.0,deg,(n-2)/2));
    return factor*PolynomialBases.gegenbauer(x,deg,(n-2)/2);
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
    z,wz = FastGaussQuadrature.gaussjacobi(4*m+deg+1,(n-3)/2,(n-3)/2);
    indices = collect(Base.product(0:1:m, 0:1:m));
    χ = ((x,i,j) -> (reproducinggegenbauer(x,n,i)*reproducinggegenbauer(x,n,j)*h(x)) );
    return map(y->( dot(wz,χ.(z,y[1],y[2]))),indices);
end

function Amatrix(n,m,k,l)
    if m < 0
        throw(DomainError(m,"order must be non-negative"));
    end
    if n < 1
        throw(DomainError(n,"dimension must be positive"));
    end
    if k > m
        throw(DomainError(k,"degree must be lesser than or equal to m"));
    end 
    if (l < 0) 
        throw(DomainError(l,"index out of bounds"))
    end
    @polyvar s
    h = reproducinggegenbauer(s,n,l);
    #TODO Fix - Horrible
    if l==0
        h = h*s^0;
    end 
    return toeplizmatrix(n,m,k,h);
end


function ffeigenvalues(n,m,k)
    invA0 = LinearAlgebra.inv(Amatrix(n,m,k,0));
    T = (1/k)*( I + sum( [ ((PolynomialBases.gegenbauer(1.0,2*j,(n-2)/2)^2 * normalizationfactor(n,0) )/((binomial(n+2*j-1,2*j)-binomial(n+2*j-3,2*j-2))^2 * normalizationfactor(n,2*j))) * invA0 * Amatrix(n,m,k,2*j) for j in 1:1:k ] ));
    V = eigen(T, sortby=abs);
    e = V.vectors[:,m+1]/norm(V.vectors[:,m+1]);
    eta = (sqrt(normalizationfactor(n,0))/surfaceareasphere(n)) * LinearAlgebra.inv(sqrt(Amatrix(n,m,k,0))) * e
    values = [];
    for i in 0:1:k
        lambda = dot(eta,Amatrix(n,m,k,2*i)*eta)*(PolynomialBases.gegenbauer(1.0,2*i,(n-2)/2)^2 * SphericalQuadrature.surfaceareasphere(n)^2)/((binomial(n+2*i-1,2*i)-binomial(n+2*i-3,2*i-2))^2 * normalizationfactor(n,2*i));
        push!(values, lambda);
    end
    return values;
end

function ffevaluation(m,p,vars)
    desc = harmonicdecomposition(p,vars)
    n = length(vars)
    k = floor(Int,maxdegree(p)/2)
    if 2*k != maxdegree(p)
        throw(DomainError(p,"polynomial p must have even degree"))
    end
    if 2*k != mindegree(p)
        throw(DomainError(p,"polynomial p must be homogeneous"))
    end
    eigvals = ffeigenvalues(n,m,2*k)
    norm = sum(vars .* vars)
    return sum( eigvals[floor(Int,u[1]/2)+1]*norm^(floor(Int,u[1]/2))*u[2] for u in desc )
end

function ffinverseeval(m,p,vars)
    desc = harmonicdecomposition(p,vars)
    n = length(vars)
    k = floor(Int,maxdegree(p)/2)
    if 2*k != maxdegree(p)
        throw(DomainError(p,"polynomial p be of even degree"))
    end
    if 2*k != mindegree(p)
        throw(DomainError(p,"polynomial p must be homogeneous"))
    end
    eigvals = ffeigenvalues(n,m,k)
    norm = sum(vars .* vars)
    return sum( (1/eigvals[floor(Int,u[1]/2)+1])*norm^(floor(Int,u[1]/2))*u[2] for u in desc )
end

