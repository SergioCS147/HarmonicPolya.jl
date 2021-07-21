module GammaOperator;

using LinearAlgebra, HarmonicBasis

function eigenvalues(n,m,k)
    if m < 0
        throw(DomainError(m,"order must be non-negative"));
    end
    if n < 1
        throw(DomainError(n,"dimension must be positive"));
    end
    if k > m
        throw(DomainError(k,"degree must be lesser than m"));
    end 
    # for i=0, the eigenvalue is always 0
    values = [1.0];
    # i>1
    for i in 1:1:k
        # might require using log to calculate in case the values get too extreme
        lambda = (2*(m-i)/(2*m+2*i+n))*values[i];
        if m==i
            lambda = factorial(m)/prod([(j+(2*m+n)/2) for j in 1:1:m]);
        end
        push!(values,lambda);
    end
    return values
end

function evaluation(m,p,x)
    desc = HarmonicBasis.harmonicdecomposition(p,x)
    n = length(x)
    k = maxdegree(p)
    if k != mindegree(p)
        throw(DomainError(p,"polynomial p must be homogeneous"))
    end
    eigvals = eigenvalues(n,m,k)
    return sum( eigvals[floor(Int,u[1]/2)+1]*norm^(floor(Int,u[1]/2))*u[2] for u in desc )
end

function inverseeval(m,p,x)
    desc = HarmonicBasis.harmonicdecomposition(p,x)
    n = length(x)
    k = maxdegree(p)
    if k != mindegree(p)
        throw(DomainError(p,"polynomial p must be homogeneous"))
    end
    eigvals = eigenvalues(n,m,k)
    return sum( (1/eigvals[floor(Int,u[1]/2)+1])*norm^(floor(Int,u[1]/2))*u[2] for u in desc )
end

end