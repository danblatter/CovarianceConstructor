# This is the main function that builds the covariance from a set of model mesh element centroids
# and a correlation kernel

using LinearAlgebra, SparseArrays
include("kernels.jl")
include("makeTearArray.jl")

function buildCovariance(C,kernel,l,tear)

T = makeTearArray(C,tear)

if cmp(kernel,"GaspariCohn") == 0
    B = buildGaspariCohn(C,l)
elseif cmp(kernel,"MattiSpecial_GC") == 0
    B = buildMattiSpecial(C,l,T)
elseif cmp(kernel,"squaredExponential") == 0
    B = buildSquaredExponential(C,l)
elseif cmp(kernel,"Exponential") == 0
    B = buildExponential(C,l)
else
    println("we only support Gaspari-Cohn and MattiSpecial_GC at the moment. Sorry!")
end

# check to see if the covariance matrix is symmetric and positive definite
if issymmetric(B)
    println("Test passed: Covariance matrix is symmetric")
else
    println("Warning! Covariance matrix is not symmetric!!!")
end
if isposdef(B)
    println("Test passed: Covariance matrix is positive definite")
else
    println("Warning! Covariance matrix is not positive definite!!!")
end

return B, T#, M, N, V

end

