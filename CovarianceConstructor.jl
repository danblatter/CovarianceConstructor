# The beginning of a new code project...

using DelimitedFiles, LinearAlgebra
include("buildCovariance.jl")

C = readdlm("resistivityMeshCentroids_Gemini.txt")  # model mesh element (model parameter) locations 

kernel = "GaspariCohn"  # correlation kernel

l = 500     # correlation length (m)

B,I,J,V = buildCovariance(C,kernel,l)