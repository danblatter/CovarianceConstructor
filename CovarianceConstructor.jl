# The beginning of a new code project...

using DelimitedFiles, JLD, SuiteSparse, PyPlot
include("buildCovariance.jl")
include("getCovSquareRoot.jl")

C = readdlm("resistivityMeshCentroids_Gemini.txt")  # model mesh element (model parameter) locations
n = size(C,1) 

kernel = "GaspariCohn"  # correlation kernel

l = 1000     # correlation length (m)

B = buildCovariance(C,kernel,l)

save("covB.jld","B",B)

L = getCovSquareRoot(B,"false")

# draw a random sample from N(0,B)
ξ = randn(n)
θ = L*ξ
# θ_noperm = L_noperm*ξ


# plot this sample

figure(1)
scatter(C[:,1],C[:,2],c=θ,s=2)
plt.xlim([-2e4, 2e4])
plt.ylim([0, 2e4])
plt.gca().invert_yaxis()

#= figure(2)
scatter(C[:,1],C[:,2],c=θ_noperm,s=2)
plt.xlim([-2e4, 2e4])
plt.ylim([0, 2e4])
plt.gca().invert_yaxis() =#

figure(3)
scatter(C[:,1],C[:,2],c=ξ,s=2)
plt.xlim([-2e4, 2e4])
plt.ylim([0, 2e4])
plt.gca().invert_yaxis()





