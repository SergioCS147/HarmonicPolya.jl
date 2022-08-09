flatten = A -> collect(Iterators.flatten(A));

function volumesphere(n)
    V = 1;
    if n==1
        V = 2;
    elseif n>1
        V = (2*π/n)*volumesphere(n-2);
    elseif n < 0
        throw(DomainError(n,"dimension argument must be non-negative"));
    end
    return V
end

function surfaceareasphere(n) 
    return (n+1)*volumesphere(n+1)
end

function sphericalquadrature(n,deg)
    x,wx = Float64[1,-1],Float64[1,1];
    if n > 1
        z,wz = gaussjacobi(deg,(n-3)/2,(n-3)/2);
        y,wy = sphericalquadrature(n-1,deg);
        xA = map(x->flatten([x[1],sqrt(1-x[1]^2).*x[2]]),Base.product(z,y))
        x = xA[:];
        wxA = map(prod,Base.product(wz,wy))
        wx = wxA[:];
    elseif n < 1
        throw(DomainError(n,"dimension argument must be positive"));
    end
    return x,wx
end