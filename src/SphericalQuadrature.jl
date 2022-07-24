flatten = A -> collect(Iterators.flatten(A));

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