using LinearAlgebra, DynamicPolynomials, Plots
include("../SphereMinimization.jl")
using .SphereMinimization

@polyvar x[1:5]
@polyvar y[1:3]
motzkin = y[1]^2*y[2]^4+ y[1]^4*y[2]^2+y[3]^6-3*y[1]^2*y[2]^2*y[3]^2
robinson = x[1]^2 * x[2]^2 + x[1]^2 * x[3]^2 + x[2]^2 * x[3]^2 + x[4]^4 - 4*x[1]*x[2]*x[3]*x[4]
vals = [SphereMinimization.solvehomogeneousgamma(motzkin,y,6)]
#vals1 = [SphereMinimization.solvehomogeneousgamma(robinson,x,6)]
vals2 = [SphereMinimization.solvehomogeneousfawzi(motzkin,y,6)]
#vals3 = [SphereMinimization.solvehomogeneousfawzi(robinson,x,6)]
for i in 4:1:50
    #push!(vals1,SphereMinimization.solvehomogeneousgamma(robinson,x,2*i));
    push!(vals,SphereMinimization.solvehomogeneousgamma(motzkin,y,2*i));
    #push!(vals3,SphereMinimization.solvehomogeneousfawzi(robinson,x,2*i));
    push!(vals2,SphereMinimization.solvehomogeneousfawzi(motzkin,y,2*i));
    if i %5 == 0
        println(i)
    end
end
plot([2*i for i in 3:1:50], [vals,vals2], label=["Motzkin - Blekherman"  "Motzkin - Fawzi" ] )
plot!(xaxis=("2s"))
plot!(yaxis=("min"))
#savefig("polys.png")
#plot!(yaxis=("min",:log10))
#savefig("polyslog.png")