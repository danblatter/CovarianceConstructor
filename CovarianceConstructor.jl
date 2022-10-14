# The beginning of a new code project...

using DelimitedFiles
include("buildCovariance.jl")

C = readdlm("resistivityMeshCentroids_Gemini.txt")

kernel = "GaspariCohn"

B = buildCovariance(C,kernel)