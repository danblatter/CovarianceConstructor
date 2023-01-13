
using PyPlot, LinearAlgebra, SparseArrays
include("kernels.jl")
include("getCovSquareRoot.jl")

myGrid = collect(0:100)
nx = length(myGrid)
X = zeros(nx,1) .+ myGrid
# l = 10 .* ones(nx)
l = [10.0]
stds = [0.3; 1; 3]

for j=1:3
    # s = stds[j] .* ones(nx)
    # B = buildMattiSpecial(X,l,s)
    B = buildGaspariCohn(X,l)
    B = stds[j] .* B
    figure(1)
    plot(B[1,:])
    L = getCovSquareRoot(B,"false")
    # draw a sample
    ξ = randn(nx)
    θ = L*ξ
    figure(2)
    plot(X,θ)
end