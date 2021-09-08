using LinearAlgebra, DynamicPolynomials, Plots
include("../SphereMinimization.jl")
using .SphereMinimization

@polyvar y[1:2]
@polyvar x[1:2]
y2 = y[1]^2
mx2 = -x[2]^2
vals = [SphereMinimization.solvehomogeneous(y2,y,6)]
vals1 = [SphereMinimization.solvehomogeneous(mx2,x,6)]
for i in 4:1:25
    push!(vals1,SphereMinimization.solvehomogeneous(mx2,x,2*i));
    push!(vals,SphereMinimization.solvehomogeneous(y2,y,2*i));
    if i %5 == 0
        println(i)
    end
end
plot([2*i for i in 3:1:25], [vals,vals1], label=["x^2" "-y^2"] )
plot!(xaxis=("2s"))
plot!(yaxis=("min"))
#savefig("graphs/monos.png")
#plot!(yaxis=("min",:log10))
#savefig("graphs/monoslog.png")