# This function will be for assigning a geologic unit to each model parameter on the basis
# of a series of horizons. The horizons are simply defined as a series of nodes (x,z) that span 
# the model from left to right in x, such that x1 is the smallest centroid x-value and xfinal is the 
# largest centroid x-value. The horizon functions are assumed piece-wise linear between nodes.
#
# For each model parameter, this function starts from the shallowest horizon (they are not allowed 
# to intersect) and proceeds to the deepest horizon until it finds a horizon that is deeper than the 
# z-coordinate of the model parameter in question. It then assigns that geologic unit's number (1:n+1,
# where n is number of horizons) to the model parameter
#

include("linearInterpolate.jl")

function assignGeologicUnits(H,C)

    # H is an array of horizons. Each horizon consists of two columns of this matrix; column i = x, 
    # column j = z, where i is odd and j is even. The deepest horizon is listed first, followed by
    # the next shallowest, until the shallowest horizon (listed last)
    # C is the array of (x,z) coordinates of model parameter centroids. First column is x-coordinates,
    # second column is z-coordinates
    # GU is the output: a vector of the same length as the model with an integer at each value
    # corresponding to the geologic unit that model parameter resides within

    # make sure there are an even number of columns in H, and determine how many horizons
    if mod(size(H,2),2) != 0
        println("There's an odd number of columns in your horizon file")
        println("Columns in this file should be in pairs (x first, then z)")
    else
        nh = Int64(size(H,2))
    end    

    # the "Geologic Unit" vector; holds the number of the geologic unit each model parameter resides
    # within
    GU = Int64.(zeros(size(C,1)))

    for ic=1:size(C,1)
        if mod(ic,100) == 0
            println("computing $ic of $(size(C,1))")
        end
        foundGeologicUnit = false
        ih = nh - 1
        while !foundGeologicUnit
        # for ih=1:2:nh-1
            # pull out this horizon only
            h = H[:,ih:ih+1]
            # append the min and max values of x to the horizon (extend horizons to the model edges) 
            firstnode = [minimum(C[:,1])-1 h[1,2]]; lastnode = [maximum(C[:,1])+1 h[end,2]];
            h = [minimum(C[:,1]) h[1,2]; h; maximum(C[:,1]) h[end,2]]
            # find the nearest two horizon nodes
            a = findmin(abs.(h[:,1] .- C[ic,1]))
            if C[ic,1] < h[a[2],1]
                ind1 = a[2] - 1; ind2 = a[2];
            elseif C[ic,1] > h[a[2],1]
                ind1 = a[2]; ind2 = a[2] + 1;
            elseif C[ic,1] == h[a[2],1]
                ind1 = a[2]; ind2 = a[2];
            end
            # setup the linear interpolation
            x1 = h[ind1,1]; x2 = h[ind2,1]; z1 = h[ind1,2]; z2 = h[ind2,2]; x = C[ic,1]; z = C[ic,2]
            zinterp = linearInterpolate(z1,z2,x1,x2,x)
            # determine if this model parameter is in this geologic unit
            if z <= zinterp
                GU[ic] = Int64((ih+1)/2 + 1)
                foundGeologicUnit = true
                # break
            elseif z > zinterp && ih == 1
                # we're above the last horizon
                GU[ic] = Int64(1)
                foundGeologicUnit = true
            else
                # move to the next horizon
                ih = ih - 2
                # debugging
                if ih < 1
                    println("z=$z; zinterp=$zinterp; ic=$ic")
                end
            end
        end
    end

    return GU

end