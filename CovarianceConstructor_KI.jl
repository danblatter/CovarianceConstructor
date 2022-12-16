# use to create covariance matrix for the King Island model

using DelimitedFiles, JLD, SuiteSparse, PyPlot
include("buildCovariance.jl")
include("getCovSquareRoot.jl")
r
C = readdlm("KI_60centroids_fullyMeshed.txt")  # model mesh element (model parameter) locations
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