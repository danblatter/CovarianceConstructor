# This is the main function that builds the covariance from a set of model mesh element centroids
# and a correlation kernel

using LinearAlgebra, SparseArrays
include("kernels.jl")

function buildCovariance(C,kernel,l)

n = size(C,1)

M = zeros(Int64,n*n,1)
N = zeros(Int64,n*n,1)
V = zeros(n*n,1)

if cmp(kernel,"GaspariCohn") == 0
    println("lets compute using the Gaspari-Cohn kernel!")
    k = 1
    for i=1:n
        if mod(i,100) == 0
            println("$i of $n")
        end
        for j=1:n
            r = norm(C[i,:] - C[j,:])
            c = GaspariCohn(r,l)
            if c > 0
                M[k] = Int(i)
                N[k] = Int(j)
                V[k] = c
                k = k + 1
            end
        end
    end
else
    println("we only support Gaspari-Cohn at the moment. Sorry!")
end

M = M[1:k-1]
N = N[1:k-1]
V = V[1:k-1]

B = sparse(M,N,V,n,n)

sparseratio = size(M,1)/(n^2)
println("This sparse covariance takes $(100*sparseratio)% of the memory of the dense matrix")

return B#, M, N, V

end

