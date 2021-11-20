using LinearAlgebra, DynamicPolynomials, Plots
include("../SphereMinimization.jl")
using .SphereMinimization

@polyvar y[1:3]
@polyvar x[1:3]
y2 = y[1]^2
mx2 = -x[2]^2
#try
    vals = [SphereMinimization.solvehomogeneousgamma(y2,y,6)]
    #vals1 = [SphereMinimization.solvehomogeneousgamma(mx2,x,6)]
    vals2 = [SphereMinimization.solvehomogeneousfawzi(y2,y,6)]
    #vals3 = [SphereMinimization.solvehomogeneousfawzi(mx2,x,6)]
#catch
    a=0
#end
for i in 4:1:75
    #try
        #push!(vals1,SphereMinimization.solvehomogeneousgamma(mx2,x,2*i));
        push!(vals,SphereMinimization.solvehomogeneousgamma(y2,y,2*i));
        #push!(vals3,SphereMinimization.solvehomogeneousfawzi(mx2,x,2*i));
        push!(vals2,SphereMinimization.solvehomogeneousfawzi(y2,y,2*i));
    #catch 
    #end
    if i %5 == 0
        println(i)
    end
end
plot([2*i for i in 3:1:75], [vals,vals2], label=["x^6, gamma" "x^6, fawzi"] )
plot!(xaxis=("2s"))
plot!(yaxis=("min"))
#savefig("graphs/monos.png")
#plot!(yaxis=("min",:log10))
#savefig("graphs/monoslog.png")