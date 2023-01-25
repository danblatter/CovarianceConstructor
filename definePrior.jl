
module definePrior

using Statistics

export assignGeologicUnits_modelParameters
export assignGeologicUnits_wellLog
export assignMeanStdCorrlen

# the goal of this function is to produce a mean and standard deviation model as well as a correlation length
# model on the basis of a well log. The basic idea is to take the well log intervals (broken up according to
# the geologic units/seismic horizons) and propagate them away from the well log into their respective geologic
# units. For each model parameter, this amounts to interpolating to a shifted and stretched version of the 
# associated well log interval mean and standard deviation

function assignMeanStdCorrlen(WL,GU,C,H)

    # WL is the well log, broken into geologic units
    #   -- abstract array with nWL elements equal to num of geologic units (1 + num of seismic horizons)
    # GU is the geologic unit designation of each model parameter
    # C contains the (x,z) locations of all the model paramters
    # H contains the seismic horizons that define the geologic units
    #   -- there are 2 columns for each horizon: x, then z
    #   -- each horizon has a unique number of nodes
    #   -- the first horizon (columns 1:2) is the surface elevation
    # z0 is the surface elevation

    meanRho = zeros(size(GU))
    stdRho = zeros(size(GU))

    nh = length(H)

    for ih=2:nh #ih=3:2:nh+2
        println("geologic unit $(ih-1) of $(nh-1)")
        # pull out the horizons bounding this geologic unit
        if ih > nh
            h_l = [0 maximum(C[:,2])]
        else
            h_l = H[ih]
        end
        h_u = H[ih-1]
        # append the min and max values of x to the horizon (extend horizons to the model edges)
        y1 = [minimum(C[:,1]) h_l[1,2]]; y2 = [maximum(C[:,1]) h_l[end,2]]
        h_l = [y1; h_l; y2]
        y1 = [minimum(C[:,1]) h_u[1,2]]; y2 = [maximum(C[:,1]) h_u[end,2]]
        h_u = [y1; h_u; y2]
        # the geologic unit these horizons correspond to
        iGU = ih-1
        # now perform the extrapolation/interpolation for every model parameter within this geologic unit
        for im in eachindex(GU)
            if GU[im] != iGU
                # skip this model parameter, as it is not within this geologic unit
                continue
            end
            # 1. Stretch well log to match geologic interval at this model parameter's location
            # find geologic interval thickness at this model parameter's location
            z_u = findHorizonDepth(C[im,1],h_u)
            z_l = findHorizonDepth(C[im,1],h_l)
            Δz_model = z_l - z_u        # z is positive down; let's keep Δz positive, too
            # well geologic interval thickness at the well location
            Δz_well = WL[iGU][end,2] - WL[iGU][1,2]
            # stretch factor
            β = Δz_model/Δz_well
            # stretch well log
            Z_well = WL[iGU][:,2] .* β
            # 2. Shift well log to match geologic interval at this model parameter's location
            # shift amount
            z_shift = Z_well[1] - z_u
            # perform shift
            Z_well = Z_well .- z_shift
            # 3. Interpolate to shifted well log (z,ρ)
            # find the two nearest well log values to this model parameter (in depth)
            ind1, ind2 = findNearestNodes(Z_well,C[im,2])
            # linear interpolation to find mean ρ at this location
            z1 = Z_well[ind1]; z2 = Z_well[ind2]; ρ1 = WL[iGU][ind1,3]; ρ2 = WL[iGU][ind2,3];
            meanRho[im] = linearInterpolate(ρ1,ρ2,z1,z2,C[im,2])
            # linear interpolation to find std ρ at this location
            ρ1 = WL[iGU][ind1,4]; ρ2 = WL[iGU][ind2,4];
            stdRho[im] = linearInterpolate(ρ1,ρ2,z1,z2,C[im,2])
            # # 4. Assign correlation length
            # if GU[im] == 1
            #     corrLen[im] = 150
            # elseif GU[im] == maximum(GU)
            #     corrLen[im] = 150
            # else
            #     corrLen[im] = 150
            # end
        end
    end

    return meanRho, stdRho

end

# This function chops a well log up into chunks based on geologic unit (horizon depths) 
#
# Input: 
# Well log array that has 4 columns
#   1. x (m) (the well could be deviated)
#   2. z (m)
#   3. mean log10 resistivity (ohm-m) (mean of log-resistivity)
#   4. log10 standard deviation (ohm-m) (standard deviation of log-resistivity)
#
# Horizon array that has an even number of columns (same as used in assignGeologicUnits)
#
# z0, the surface elevation at the well location (m)
#
# zmax, the bottom of the region of interest (m); the depth to which you want well log values to be 
# interpolated to. Below this, the model is assumed to be constant
#
# Output: an abstract array of length equal to the number of geologic units. Each element of this
# array is itself an array with the same 4 columns as the input, but where the number of rows is 
# equal to the number of well log measurements lying within that geologic unit. In addition to the
# well log measurements, the first row of each 4-column sub-array will be the top depth of the geologic
# unit and the last row will be the bottom depth of the unit. The x, mean(rho), and std(rho) values
#  of the top and bottom rows will be copies of the second and next-to-last rows, respectively.

function assignGeologicUnits_wellLog(wellLog,H,zmax)

    nh = length(H)

    wellLogUnits = Array{Float64}[]     # construct an as-yet empty vector of arrays
                                        # each array in this vector will be a well log unit

    global zinterp = 0  # placeholder value, just declaring it here so it persists beyond while loop
    ih = 2
    while ih <= nh   
        println("ih = $ih")
        # pull out this horizon only
        h = H[ih]
        # x-location of the well (assume it's not deviated)
        x_well = mean(wellLog[:,1])
        # find the nearest two horizon nodes to the well 
        ind1,ind2 = findNearestNodes(h[:,1],x_well)
        # find depth to this horizon at location of well (assume it's not deviated for now)
        x1 = h[ind1,1]; x2 = h[ind2,1]; z1 = h[ind1,2]; z2 = h[ind2,2]; x = x_well;
        if ih == 2 
            # shallowest horizon; zprev is the top of the model
            h_ = H[1]
            # find surface elevation at location of well (assume it's not deviated for now)
            ind1,ind2 = findNearestNodes(h_[:,1],x_well)
            z1_ = h_[ind1,2]; z2_ = h_[ind2,2];
            zprev = linearInterpolate(z1_,z2_,x1,x2,x)
            global zinterp = linearInterpolate(z1,z2,x1,x2,x)
            println("ih=$ih, z_upper=$zprev, z_lower=$zinterp")
        else
            zprev =  zinterp    # depth of previous horizon
            global zinterp = linearInterpolate(z1,z2,x1,x2,x)
            println("ih=$ih, z_upper=$zprev, z_lower=$zinterp")
        end
        # find all the well log values in this geologic unit
        p = findall(t -> t .< zinterp, wellLog[:,2])
        q = findall(t -> t .> zprev, wellLog[:,2])
        inds = intersect(p,q)
        println("Looking between $zprev and $zinterp, we have $(length(inds)) well log values")
        tmp = wellLog[inds,:]
        # tack an upper and lower row to represent the bounds of this geologic unit
        thisUnit = [tmp[1,1] zprev tmp[1,3:4]'; tmp; tmp[end,1] zinterp tmp[end,3:4]']
        # add this well log section to our vector of arrays
        push!(wellLogUnits,thisUnit)

        ih = ih + 1

    end

    return wellLogUnits

end

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

function assignGeologicUnits_modelParameters(H,C)

    # H is an array of arrays, each of which has two columns: column 1 is x-coordinates in meters
    # while column 2 is z-coordinates, positive down, in meters.

    # H is an array of horizons. Each horizon consists of two columns of this matrix; column i = x, 
    # column j = z, where i is odd and j is even. The shallowest horizon is listed first, followed by
    # the next deepest, until the deepest horizon (listed last)
    # C is the array of (x,z) coordinates of model parameter centroids. First column is x-coordinates,
    # second column is z-coordinates
    # GU is the output: a vector of the same length as the model with an integer at each value
    # corresponding to the geologic unit that model parameter resides within

    nh = length(H)
    println("we have $nh geologic surfaces, including top and bottom of model")

    # the "Geologic Unit" vector; holds the number of the geologic unit each model parameter resides
    # within
    GU = Int64.(zeros(size(C,1)))

    for ic=1:size(C,1)
        if mod(ic,100) == 0
            println("computing $ic of $(size(C,1))")
        end
        foundGeologicUnit = false 
        ih = 2      # first surface must be the upper boundary of the model
        while !foundGeologicUnit
            # # pull out this horizon only
            h = H[ih]
            # append the min and max values of x to the horizon (extend horizons to the model edges) 
            h = [minimum(C[:,1]) h[1,2]; h; maximum(C[:,1]) h[end,2]]
            # find the nearest two horizon nodes
            ind1,ind2 = findNearestNodes(h[:,1],C[ic,1])
            # setup the linear interpolation
            x1 = h[ind1,1]; x2 = h[ind2,1]; z1 = h[ind1,2]; z2 = h[ind2,2]; x = C[ic,1]; z = C[ic,2]
            # this is the horizon depth at this model parameter's x-location
            zinterp = linearInterpolate(z1,z2,x1,x2,x) 
            # determine if this model parameter is in this geologic unit
            if z <= zinterp
                GU[ic] = ih - 1
                foundGeologicUnit = true
            elseif z > zinterp && ih == nh
                # we're below the last horizon
                GU[ic] = nh
                foundGeologicUnit = true
            else 
                # move to the next horizon
                ih = ih + 1
            end
        end
    end

    return GU

end

function findHorizonDepth(x,h)
    # find the depth to a horizon, for a given model parameter x-location
    # x is the model parameter x-location (m)
    # h is the horizon: rows are the (x,z nodes); first column of h is x, second column is z (m)

    # find the nearest two horizon nodes
    a = findmin(abs.(h[:,1] .- x))
    if x < h[a[2],1]
        ind1 = a[2] - 1; ind2 = a[2];
    elseif x > h[a[2],1]
        ind1 = a[2]; ind2 = a[2] + 1;
    elseif x == h[a[2],1]
        ind1 = a[2]; ind2 = a[2];
    end
    # setup the linear interpolation
    x1 = h[ind1,1]; x2 = h[ind2,1]; z1 = h[ind1,2]; z2 = h[ind2,2];
    # this is the horizon depth at this model parameter's x-location
    zinterp = linearInterpolate(z1,z2,x1,x2,x) 
    return zinterp

end

function findNearestNodes(X,x)
    # This function finds the nearest values in an array X to a point x
    a = findmin(abs.(X .- x))
    if x < X[a[2]]
        ind1 = a[2] - 1; ind2 = a[2];
    elseif x > X[a[2]]
        ind1 = a[2]; ind2 = a[2] + 1;
    elseif x == X[a[2]]
        ind1 = a[2]; ind2 = a[2];
    end
    return ind1, ind2
end

function linearInterpolate(y0,y1,x0,x1,x)

    if x0 == x1
        y = y0
    else
        y = y0 + (x-x0)*(y1-y0)/(x1-x0)
    end

    return y

end

end