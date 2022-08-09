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

@testset "Quadrature Tests" begin
    @testset "Quadrature Volume on Dimension" begin 
        for n in 1:10
            z,w = HarmonicPolya.sphericalquadrature(n+1,2)
            V = sum(w)
            V1 = surfaceareasphere(n)
            @test V ≈ V1 atol=0.001
        end
    end 

    @testset "Quadrature Volume on Degree" begin 
        V2 = surfaceareasphere(2)

        function testquadrature()
            for d in 2:10
                z1,w1 = HarmonicPolya.sphericalquadrature(3,2*d)
                V3 = sum(w1)
                @test V2 ≈ V3 atol=0.001
            end
        end
    end

    #@testset "Quadrature Test on Monomials" begin
    #    a=0
    #end
end