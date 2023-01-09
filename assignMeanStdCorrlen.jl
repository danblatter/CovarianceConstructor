
# module assignPriorInformation

include("linearInterpolate.jl")

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

    for ih=3:2:nh+2
        println("computing horizon $(floor((ih+1)/2)-1) of $(floor(nh/2)-1)")
        # pull out the horizons bounding this geologic unit
        if ih > nh
            h_l = [0 maximum(C[:,2])]
        else
            inds = findall(x -> typeof(x) == Float64, H[:,ih])
            h_l = H[inds,ih:ih+1]
        end
        inds = findall(x -> typeof(x) == Float64, H[:,ih-2])
        h_u = H[inds,ih-2:ih-1]
        # append the min and max values of x to the horizon (extend horizons to the model edges)
        y1 = [minimum(C[:,1]) h_l[1,2]]; y2 = [maximum(C[:,1]) h_l[end,2]]
        h_l = [y1; h_l; y2]
        y1 = [minimum(C[:,1]) h_u[1,2]]; y2 = [maximum(C[:,1]) h_u[end,2]]
        h_u = [y1; h_u; y2]
        # the geologic unit these horizons correspond to
        iGU = Int64((ih-1)/2)   # GU 1 corresponds to horizon 1, which is H[:,3:4] because H[:,1:2] is the surface
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
            Δz_well = wellLogUnits[iGU][end,2] - wellLogUnits[iGU][1,2]
            # stretch factor
            β = Δz_model/Δz_well
            # stretch well log
            Z_well = wellLogUnits[iGU][:,2] .* β
            # 2. Shift well log to match geologic interval at this model parameter's location
            # shift amount
            z_shift = Z_well[1] - z_u
            # perform shift
            Z_well = Z_well .- z_shift
            # # at this point, z_u and z_l should match Z_well[1] and Z_well[end]
            # println("beta = $β")
            # println("transformed well interval bounds: ($(Z_well[1]), $(Z_well[end])); ")
            # println("model location interval bounds: ($z_u, $z_l)")
            # println(" ")
            # 3. Interpolate to shifted well log (z,ρ)
            # find the two nearest well log values to this model parameter (in depth)
            ind1, ind2 = findNearestNodes(Z_well,C[im,2])
            # linear interpolation to find mean ρ at this location
            z1 = Z_well[ind1]; z2 = Z_well[ind2]; ρ1 = wellLogUnits[iGU][ind1,3]; ρ2 = wellLogUnits[iGU][ind2,3];
            meanRho[im] = linearInterpolate(ρ1,ρ2,z1,z2,C[im,2])
            # linear interpolation to find std ρ at this location
            ρ1 = wellLogUnits[iGU][ind1,4]; ρ2 = wellLogUnits[iGU][ind2,4];
            stdRho[im] = linearInterpolate(ρ1,ρ2,z1,z2,C[im,2])
            # 4. Assign correlation length to be interval thickness at model parameter location
            if GU[im] == 1
                corrLen[im] = 500
            elseif GU[im] == maximum(GU)
                corrLen[im] = 1500
            else
                corrLen[im] = Δz_model
            end
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

# end