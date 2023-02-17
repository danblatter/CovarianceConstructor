# This is the main function that builds the covariance from a set of model mesh element centroids
# and a correlation kernel
# The input arguments are 
#   C: an nm x 2 array of model parameter centroids (x,z), where nm is # of model parameters
#   kernel: a string specifying the correlation kernel to use
#   l: the correlation length; could be a scalar, a function, or an nm x 2 array
#   s: the model standard deviation; can be a scalar or an nm x 2 array

using LinearAlgebra, SparseArrays
include("kernels.jl")

function buildCovariance(C,kernel,l,s)

if cmp(kernel,"GaspariCohn") == 0
    B = buildGaspariCohn(C,l)
elseif cmp(kernel,"MattiSpecial_GC") == 0
    B = buildMattiSpecial(C,l,s)
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

return B #, M, N, V

end

