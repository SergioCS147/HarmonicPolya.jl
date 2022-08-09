@testset "Minimization Benchmarks" begin
    @polyvar x[1:4]
    robinson = x[1]^2*(x[1]-x[4])^2 + x[2]^2*(x[2]-x[4])^2 + x[3]^2*(x[3]-x[4])^2 + 2*x[1]*x[2]*x[3]*(x[1]+x[2]+x[3]-2*x[4])
    @polyvar y[1:3]
    motzkin = y[1]^2*y[2]^4+ y[1]^4*y[2]^2+y[3]^6-3*y[1]^2*y[2]^2*y[3]^2
    @test HarmonicPolya.lowerboundsquares(robinson,x,50) ≈ 0 atol=0.5
    @test HarmonicPolya.lowerboundfawzi(robinson,x,50) ≈ 0 atol=0.025
    @test HarmonicPolya.lowerboundsquares(motzkin,y,50) ≈ 0 atol=0.5
    @test HarmonicPolya.lowerboundfawzi(motzkin,y,50) ≈ 0 atol=0.025
end