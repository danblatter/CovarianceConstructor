# This function will be for assigning a geologic unit to each model parameter on the basis
# of a series of horizons. The horizons are simply defined as a series of nodes (x,z) that span 
# the model from left to right in x, such that x1 is the smallest centroid x-value and xfinal is the 
# largest centroid x-value. The horizon functions are assumed piece-wise linear between nodes.
#
# For each model parameter, this function starts from the shallowest horizon (they are not allowed 
# to intersect) and proceeds to the deepest horizon until it finds a horizon that is deeper than the 
# z-coordinate of the model parameter in question. It then assigns that geologic unit's number (1:n+1,
# where n is number of horizons) to the model parameter

function assignGeologicUnits(H,m)

    # H is an array of horizons. Each horizon consists of two columns of this matrix; column i = x, 
    # column j = z, where i is odd and j is even

    if mod(size(H,2),2) != 0
        println("There's an odd number of columns in your horizon file")
        println("Columns in this file should be in pairs (x first, then z)")
    else
        nh = Int64(size(H,2)/2)
    end

end