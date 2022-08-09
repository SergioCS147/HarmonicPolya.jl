function generatemonomials(x,m,;option=0)
    mons = [];
    if option == 0
        mons = monomials(x,[1:m;], mon -> (DynamicPolynomials.degree(mon) == m));
    elseif 1 <= option <= length(x)
        mons = monomials(x,[1:m;], mon -> ((DynamicPolynomials.degree(mon) == m)&&(DynamicPolynomials.degree(mon,x[option])<=1)));
    else 
        throw(DomainError(option,"option must be non-negative and <=$(length(x))"))
    end
    return mons
end

function differentiatepow(u,var::MP.AbstractVariable;power=1)
    res = u;
    if power>1
        res = DynamicPolynomials.differentiate(differentiatepow(u,var,power=power-1),var);
    elseif power==1
        res = DynamicPolynomials.differentiate(u,var);
    elseif power < 0
        throw(DomainError(power,"power must be non-negative"));
    end
    return res;
end

function dalpha(u,alpha,vars)
    if length(alpha) != length(vars)
        throw(DomainError(alpha,"alpha must be a valid multi-index for the variables x"))
    end
    res = u;
    for i in 1:1:length(alpha)
        res = differentiatepow(res,vars[i],power=alpha[i]);
    end 
    return res;
end

function polyderivative(u,poly,vars)
    pterms = terms(poly);
    pcoefs = map(coefficient,pterms);
    pmons = map(monomial,pterms);
    fixedexp = z -> DynamicPolynomials.exponents(prod(vars.^0)*z);
    ppows = map(fixedexp,pmons);
    Dvec = map(z->dalpha(u,z,vars),ppows);
    if length(pcoefs)!= length(Dvec)
        throw(ErrorException("poly was not well defined"))
    end
    return sum(pcoefs[i]*Dvec[i] for i in 1:1:length(Dvec));
end

function laplacian(u,vars,;power=1)
    res = u;
    if power > 1
        res = polyderivative(laplacian(u,vars,power=power-1),sum(vars.*vars),vars);
    elseif power == 1
        res = polyderivative(u,sum(vars .* vars),vars);
    elseif power < 0 
        throw(DomainError(power,"power must be non-negative"));
    end 
    return res;
end

function normpartialder(hvec,var::MP.AbstractVariable,vars,;power=1)
    res = hvec;
    if power>1
        der = normpartialder(hvec,var,vars,power=power-1)
        res = (der[1]-2,der[1]*var*der[2]+sum(vars.*vars)*DynamicPolynomials.differentiate(der[2],var));
    elseif power==1
        res = (hvec[1]-2,hvec[1]*var*hvec[2]+sum(vars.*vars)*DynamicPolynomials.differentiate(hvec[2],var));
    elseif power < 0
        throw(DomainError(power,"power must be non-negative"));
    end
    return res;
    return 
end

function normdalpha(hvec,alpha,vars)
    if length(alpha) != length(vars)
        throw(DomainError(alpha,"alpha must be a valid multi-index for the variables x"))
    end
    res = hvec;
    for i in 1:1:length(alpha)
        res = normpartialder(res,vars[i],vars,power=alpha[i]);
    end 
    return res;
end

function normpolyderivative(hvec,poly,vars)
    pterms = terms(poly);
    pcoefs = map(coefficient,pterms);
    pmons = map(monomial,pterms);
    fixedexp = z -> DynamicPolynomials.exponents(prod(vars.^0)*z);
    ppows = map(fixedexp,pmons);
    Dvec = map(z->normdalpha(hvec,z,x),ppows);
    if length(pcoefs)!= length(Dvec)
        throw(ErrorException("poly was not well defined"))
    end
    return Dvec
end

function generatebasissphere(m, vars)
    n = length(vars)
    genpoly = (2-n,1);
    monos = generatemonomials(vars,m,option=1);
    fixedexp = z -> DynamicPolynomials.exponents(prod(vars.^0)*z);
    mpows = map(fixedexp,monos)
    return map(m->normdalpha(genpoly,m,vars)[2],mpows)
end

function dcoefficient(n,m,i,j)
    d = 0;
    if i == j
        d = 1.0/(2.0^j*factorial(j)*prod([(2*m + n - 2*j - 2*l) for l in 1:j]));
    elseif i < j
        d = -dcoefficient(n,m,i,j-1)/(2*(j-i)*(2*m+n-2-2*i-2*j));
    end
    return d;
end

function harmonicdecomposition(u,vars)
    norm = sum(vars .* vars);
    n =length(vars);
    m = maxdegree(u);
    if mindegree(u) != m
        throw(DomainError(u,"polynomial must be homogeneous"));
    end
    return [ (2*i,sum([dcoefficient(n,m,i,j)*norm^(j-i)*laplacian(u,vars,power=j) for j in i:floor(Int,m/2)])) for i in 0:floor(Int,m/2) ];
end

