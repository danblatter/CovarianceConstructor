# The beginning of a new code project...

using DelimitedFiles, LinearAlgebra, JLD
include("buildCovariance.jl")

C = readdlm("resistivityMeshCentroids_Gemini.txt")  # model mesh element (model parameter) locations 

kernel = "GaspariCohn"  # correlation kernel

l = 500     # correlation length (m)

B = buildCovariance(C,kernel,l)

save("covB.jld","B",B)