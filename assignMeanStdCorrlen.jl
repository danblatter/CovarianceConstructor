
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
    corrLen = zeros(size(GU))

    nh = size(H,2)

    for ih=3:2:nh
        # pull out the horizons bounding this geologic unit
        inds = findall(x -> typeof(x) == Float64, H[:,ih])
        h_l = H[inds,ih:ih+1]
        inds = findall(x -> typeof(x) == Float64, H[:,ih-2])
        h_u = H[inds,ih-2:ih-1]
        # append the min and max values of x to the horizon (extend horizons to the model edges) 
        h_l = [minimum(C[:,1]) h_l[1,2]; h_l; maximum(C[:,1]) h_l[end,2]]
        h_u = [minimum(C[:,1]) h_u[1,2]; h_u maximum(C[:,1]) h_u[end,2]]
        # the geologic unit these horizons correspond to
        iGU = Int64((ih-1)/2)   # GU 1 corresponds to horizon 1, which is H[:,3:4] because H[:,1:2] is the surface
        # now perform the extrapolation/interpolation for every model parameter within this geologic unit
        for im in eachindex(GU)
            if GU[im] != iGU
                # skip this model parameter, as it is not within this geologic unit
                continue
            end
            # find top and bottom depth of this geologic interval at this model parameter's location
            z_u = findHorizonDeph(C[im,1],h_u)
            z_l = findHorizonDeph(C[im,1],h_l)
        end
    end

    return meanRho, stdRho, corrLen

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