import HarmonicPolya

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

function testquadrature()
    for n in 1:10
        z,w = HarmonicPolya.sphericalquadrature(n,2)
        V = sum(w)
        
    end
end