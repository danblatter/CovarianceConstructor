# This is the main function that builds the covariance from a set of model mesh element centroids
# and a correlation kernel

function buildCovariance(C,kernel)

B = zeros(size(C,1),size(C,1))    
println("the size of B is $(size(B,1)) x $(size(B,2))")

if cmp(kernel,"GaspariCohn") == 0
    println("lets compute using the Gaspari-Cohn kernel!")
else
    println("we only support Gaspari-Cohn at the moment. Sorry!")
end

return B

end