# The beginning of a new code project...

using DelimitedFiles, JLD, SuiteSparse, PyPlot
include("buildCovariance.jl")
include("getCovSquareRoot.jl")

C = readdlm("resistivityMeshCentroids_Gemini.txt")  # model mesh element (model parameter) locations
n = size(C,1) 
tear = readdlm("tear.txt")
# tear = Inf                 # uncomment this line if you don't wish to add a correlation tear

kernel = "GaspariCohn"  # correlation kernel

l = 1000     # correlation length (m)

B, T = buildCovariance(C,kernel,l,tear)

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

figure(4)
scatter(C[:,1],C[:,2],c=T,s=2)
plt.xlim([-2e4, 2e4])
plt.ylim([0, 2e4])
plt.gca().invert_yaxis()
colorbar()




