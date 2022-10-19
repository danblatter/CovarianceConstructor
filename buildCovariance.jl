# This is the main function that builds the covariance from a set of model mesh element centroids
# and a correlation kernel

using LinearAlgebra
include("kernels.jl")

function buildCovariance(C,kernel,l)

n = size(C,1)

I = zeros(Int64,n*n,1)
J = zeros(Int64,n*n,1)
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
                I[k] = Int(i)
                J[k] = Int(j)
                V[k] = c
                k = k + 1
            end
        end
    end
else
    println("we only support Gaspari-Cohn at the moment. Sorry!")
end

I = I[1:k-1]
J = J[1:k-1]
V = V[1:k-1]

B = sparse(I,J,V,n,n)

sparseratio = size(I,1)/(n^2)
println("This sparse covariance takes $(100*sparseratio)% of the memory of the dense matrix")

return B, I, J, V

end

