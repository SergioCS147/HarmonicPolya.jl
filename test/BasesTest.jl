

@testset "Test Laplacian" begin
    @polyvar x[1:3]
    t1 = x[1]^2-2*x[2]^2 +x[3]^2
    t2 = x[1]^3-3*x[1]*x[3]^2

    @test HarmonicPolya.laplacian(t1, x) == 0
    @test HarmonicPolya.laplacian(t2, x) == 0
    @test HarmonicPolya.laplacian(t1+t2, x) == 0

    t3 = x[1]^2*x[2]+x[3]^3

    @test HarmonicPolya.laplacian(t3, x) == 2*x[2] + 6*x[3]
    @test HarmonicPolya.laplacian(t3+t1+t2, x) == 2*x[2] + 6*x[3]
end 

@testset "Test Harmonic Basis" begin
    @polyvar x[1:3]
    for d in 1:5
        basis = HarmonicPolya.generatebasissphere(2*d,x)
        for i in 1:length(basis)
            @test laplacian(basis[i],x) == 0
        end
    end
end