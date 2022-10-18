# The beginning of a new code project...

using DelimitedFiles, LinearAlgebra
include("buildCovariance.jl")

C = readdlm("resistivityMeshCentroids_Gemini.txt")

kernel = "GaspariCohn"

B,I,J,V = buildCovariance(C,kernel)