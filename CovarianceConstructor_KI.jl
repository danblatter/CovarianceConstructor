# use to create covariance matrix for the King Island model

using DelimitedFiles, JLD, SparseArrays, PyPlot, Statistics
include("buildCovariance.jl")
include("getCovSquareRoot.jl")
include("makeSyntheticKIwell.jl")
include("definePrior.jl")
using .definePrior

println("reading in the model centroids")
# C = readdlm("KI_60centroids_fullyMeshed.txt")  # model mesh element (model parameter) locations
# C = readdlm("KI_60centroids_justReservoir2.txt")  # model mesh element (model parameter) locations
C = readdlm("MARE2DEM_models/KIyr08_60deg_centroids_justReservoir.txt")  # model mesh element (model parameter) locations
# trueRho = readdlm("KI_60rho_fullyMeshed.txt")  # model mesh element (model parameter) locations
# trueRho = readdlm("KI_60rho_justReservoir2.txt")  # model mesh element (model parameter) locations
trueRho = readdlm("MARE2DEM_models/KIyr08_60deg_rho_justReservoir.txt")  # model mesh element (model parameter) locations
n = size(C,1) 

# kernel = "GaspariCohn"  # correlation kernel
#   Note on the GaspariCohn kernel: the correlation length can be anisotropic in any number off
#   spatial dimensions, but must be stationary. The standard deviation can be non-stationary
# kernel = "Exponential"
# kernel = "squaredExponential"
#   Note on the Matti special: it can be called with three types of arguments for "l"
#   1. a scalar; l is the same everywhere
#   2. a function; must be a function of x and z
# kernel = "MattiSpecial_GC"  # correlation kernel
# corrlen_MSgc = 150    # isotropic correlation length (m)
kernel = "GaspariCohn"  # correlation kernel
corrLen = [500.0; 25.0]     # anisotropic correlation lengths, [x ;z] (m)

# load in the well log; we'll build the prior mean and standard deviation models from this (along with 
# geologic/seismic surfaces)
println("loading well log")
wellX = 0       # injection well location (m)
wellLog = makeSyntheticKIwell(C,trueRho,wellX)

# load in the geologic horizons
println("loading geologic horizons")
s_top = readdlm("Surfaces/KI60_topSeal.txt")
# s2 = readdlm("Surfaces/KI60_topSeal_mine.txt")
# s3 = readdlm("Surfaces/KI60_topRes.txt")
s_bot = readdlm("Surfaces/KI60_baseRes.txt")
# H = Array[s_top,s2,s3,s_bot]
H = Array[s_top,s_bot]

# get the geologic unit assignment vector
println("assigning each model parameter to a geologic unit")
GU = assignGeologicUnits_modelParameters(H,C)

# get the well log units, divided up according to geologic unit
println("breaking up the well log based on geologic unit")
zmax = maximum(C[:,2])     # make sure the 'well log' extends to cover the deepest model parameter  
                           # portions of the prior below the depth of actual prior information from
                           # a well log are assumed to have the same value as the last well log value
wellLogUnits = assignGeologicUnits_wellLog(wellLog,H,zmax)

meanRho, stdRho = assignMeanStdCorrlen(wellLogUnits,GU,C,H)

if cmp(kernel,"MattiSpecial_GC") == 0
    corrLen = corrlen_MSgc .* ones(size(meanRho))
    B = buildCovariance(C,kernel,corrLen,stdRho)
else
    B = buildCovariance(C,kernel,corrLen,stdRho)
    B = (mean(stdRho))^2 .* B
end

figure(1,figsize=(6.0,2.75))
scatter(C[:,1],C[:,2],c=B[Int64(ceil(size(B,1)*rand())),:],s=5)
plt.xlim([-5e3, 7e3])
plt.ylim([1.1e3, 2.25e3])
plt.gca().invert_yaxis()
plt.ylabel("depth (m)")
plt.xlabel("horiz. position (m)")
plt.title("correlation kernel")

save("covB_KI.jld","B",B)

println("computing matrix square root...")
L = getCovSquareRoot(B,"false")

println("drawing random sample...")
# draw a random sample from N(θ_mean,B)
ξ = randn(n)
θ = L*ξ .+ meanRho

println("plotting...")
# plot this sample

figure(2,figsize=(6.5,2.75))
scatter(C[:,1],C[:,2],c=θ,s=5,cmap=ColorMap("turbo"))
plt.xlim([-5e3, 7e3])
plt.ylim([1.1e3, 2.25e3])
plt.ylabel("depth (m)")
plt.xlabel("horiz. position (m)")
plt.title("random prior ρ model")
plt.gca().invert_yaxis()
colorbar()
clim(0, 1.75)

figure(3,figsize=(6.5,2.75))
scatter(C[:,1],C[:,2],c=trueRho,s=5,cmap=ColorMap("turbo"))
plt.xlim([-5e3, 7e3])
plt.ylim([1.1e3, 2.25e3])
plt.ylabel("depth (m)")
plt.xlabel("horiz. position (m)")
plt.title("true ρ model")
plt.gca().invert_yaxis()
colorbar()
clim(0, 1.75)

figure(4,figsize=(6.5,2.75))
scatter(C[:,1],C[:,2],c=meanRho,s=5,cmap=ColorMap("turbo"))
plt.xlim([-5e3, 7e3])
plt.ylim([1.1e3, 2.25e3])
plt.ylabel("depth (m)")
plt.xlabel("horiz. position (m)")
plt.title("mean ρ model")
plt.gca().invert_yaxis()
colorbar()
clim(0, 1.75)

figure(5,figsize=(6.5,2.75))
scatter(C[:,1],C[:,2],c=GU,s=5,cmap=ColorMap("turbo"))
plt.xlim([-5e3, 7e3])
plt.ylim([1.1e3, 2.25e3])
plt.ylabel("depth (m)")
plt.xlabel("horiz. position (m)")
plt.title("geologic units (and well sections)")
plt.gca().invert_yaxis()
for ip in eachindex(wellLogUnits)
    plot(wellLogUnits[ip][:,1],wellLogUnits[ip][:,2])
end
