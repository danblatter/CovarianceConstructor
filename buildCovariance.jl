# This is the main function that builds the covariance from a set of model mesh element centroids
# and a correlation kernel

using LinearAlgebra
include("kernels.jl")

function buildCovariance(C,kernel)

#= B = zeros(size(C,1),size(C,1))    
println("the size of B is $(size(B,1)) x $(size(B,2))") =#

n = size(C,1)

I = zeros(Int64,n*n,1)
J = zeros(Int64,n*n,1)
V = zeros(n*n,1)

if cmp(kernel,"GaspariCohn") == 0
    println("lets compute using the Gaspari-Cohn kernel!")
#=     x = 0:0.01:3
    y = zeros(size(x))
    for j=1:length(x)
        y[j] = GaspariCohn(x[j],1)
    end
    writedlm("GaspariCohnTest.txt",y) =#
    k = 1
    for i=1:n
        println("$i of $n")
        for j=1:n
            r = norm(C[i,:] - C[j,:])
            c = GaspariCohn(r,1e3)
            # println("r = $r; c = $c")
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

return B, I, J, V

end

