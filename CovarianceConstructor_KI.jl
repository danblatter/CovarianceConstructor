# use to create covariance matrix for the King Island model

using DelimitedFiles, JLD, SparseArrays, PyPlot, StatsBase
include("buildCovariance.jl")
include("getCovSquareRoot.jl")
include("assignGeologicUnits_modelParameters.jl")
include("assignGeologicUnits_wellLog.jl")
include("assignMeanStdCorrlen.jl")

println("reading in the model centroids")
C = readdlm("KI_60centroids_fullyMeshed.txt")  # model mesh element (model parameter) locations
trueRho = readdlm("KI_60rho_fullyMeshed.txt")  # model mesh element (model parameter) locations
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
# kernel = "GaspariCohn"  # correlation kernel

# load in the well log
println("loading well log")
wellLog = readdlm("KIrhoWellLog.txt")

# load in the geologic horizons
println("loading geologic horizons")
H = readdlm("KI60horizons_3_wSurface.txt")

# get the geologic unit assignment vector
println("assigning each model parameter to a geologic unit")
GU = assignGeologicUnits_modelParameters(H,C)

# get the well log units, divided up according to geologic unit
println("breaking up the well log based on geologic unit")
zmax = maximum(C[:,2])     # make sure the 'well log' extends to cover the deepest model parameter  
                           # portions of the prior below the depth of actual prior information from
                           # a well log are assumed to have the same value as the last well log value
wellLogUnits = assignGeologicUnits_wellLog(wellLog,H,zmax)

meanRho, stdRho, corrLen = assignMeanStdCorrlen(wellLogUnits,GU,C,H)

if cmp(kernel,"MattiSpecial_GC") == 0
B = buildCovariance(C,kernel,corrLen,stdRho)
else
    corrLen = [1000.0; 75.0]
    B = buildCovariance(C,kernel,corrLen,stdRho)
    B = mean(stdRho) .* B
end

figure(1,figsize=(6,2.5))
scatter(C[:,1],C[:,2],c=B[4000,:],s=2)
plt.xlim([-5e3, 7e3])
plt.ylim([0, 2.25e3])
plt.gca().invert_yaxis()

save("covB_KI.jld","B",B)

println("computing matrix square root...")
L = getCovSquareRoot(B,"false")

println("drawing random sample...")
# draw a random sample from N(θ_mean,B)
ξ = randn(n)
θ = L*ξ .+ meanRho

println("plotting...")
# plot this sample

figure(2,figsize=(6,2.5))
scatter(C[:,1],C[:,2],c=θ,s=2,cmap=ColorMap("turbo"))
plt.xlim([-5e3, 7e3])
plt.ylim([0, 2.25e3])
plt.gca().invert_yaxis()
colorbar()

figure(3,figsize=(6,2.5))
scatter(C[:,1],C[:,2],c=trueRho,s=2,cmap=ColorMap("turbo"))
plt.xlim([-5e3, 7e3])
plt.ylim([0, 2.25e3])
plt.gca().invert_yaxis()
colorbar()

figure(4,figsize=(6,2.5))
scatter(C[:,1],C[:,2],c=GU,s=2,cmap=ColorMap("turbo"))
plt.xlim([-5e3, 7e3])
plt.ylim([0, 2.25e3])
plt.gca().invert_yaxis()
