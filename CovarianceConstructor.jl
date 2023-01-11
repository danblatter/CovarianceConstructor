# The beginning of a new code project...

using DelimitedFiles, JLD, SuiteSparse, PyPlot
include("buildCovariance.jl")
include("getCovSquareRoot.jl")

C = readdlm("resistivityMeshCentroids_Gemini.txt")  # model mesh element (model parameter) locations
n = size(C,1) 
# tear = readdlm("tear.txt")
tear = Inf                 # uncomment this line if you don't wish to add a correlation tear

# kernel = "GaspariCohn"  # correlation kernel
# kernel = "Exponential"
# kernel = "squaredExponential"
#   Note on the Matti special: it can be called with three types of arguments for "l"
#   1. a scalar; l is the same everywhere
#   2. a function; must be a function of x and z
kernel = "MattiSpecial_GC"  # correlation kernel

# l = 2500     # scalar correlation length (m)
l(x::Float64,z::Float64) = 250 + 0.5*z      # correlation length as a function of x and z 

s = 1.0       # leave s=1 unless using non-stationary kernel and want to specify non-stationary variance

B = buildCovariance(C,kernel,l,s)

figure(1)
scatter(C[:,1],C[:,2],c=B[4000,:],s=2)
plt.xlim([-2e4, 2e4])
plt.ylim([0, 2e4])
plt.gca().invert_yaxis()

save("covB.jld","B",B)

println("computing matrix square root...")
L = getCovSquareRoot(B,"false")

println("drawing random sample...")
# draw a random sample from N(0,B)
ξ = randn(n)
θ = L*ξ
# θ_noperm = L_noperm*ξ

println("plotting...")
# plot this sample

figure(2)
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

println("finished!")


