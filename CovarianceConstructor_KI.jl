# use to create covariance matrix for the King Island model

using DelimitedFiles, JLD, SuiteSparse, PyPlot
include("buildCovariance.jl")
include("getCovSquareRoot.jl")
include("assignGeologicUnits_modelParameters.jl")
include("assignGeologicUnits_wellLog.jl")

println("reading in the model centroids")
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

# load in the well log
println("loading well log")
wellLog = readdlm("KIrhoWellLog.txt")

# load in the geologic horizons
println("loading geologic horizons")
H = readdlm("KI60horizons_3.txt")

# get the geologic unit assignment vector
println("assigning each model parameter to a geologic unit")
GU = assignGeologicUnits_modelParameters(H,C)

# get the well log units, divided up according to geologic unit
println("breaking up the well log based on geologic unit")
wellLogUnits = assignGeologicUnits_wellLog(wellLog,H,2300,0)
