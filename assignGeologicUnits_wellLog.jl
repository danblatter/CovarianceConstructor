
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

using StatsBase, DelimitedFiles
include("linearInterpolate.jl")

function assignGeologicUnits_wellLog(wellLog,H,zmax)

   # make sure there are an even number of columns in H, and determine how many horizons
   if mod(size(H,2),2) != 0
        println("There's an odd number of columns in your horizon file")
        println("Columns in this file should be in pairs (x first, then z)")
    else
        nh = Int64(size(H,2))
    end   

    wellLogUnits = Array{Float64}[]     # construct an as-yet empty vector of arrays
                                        # each array in this vector will be a well log unit

    global zinterp = 0  # placeholder value, just declaring it here so it persists beyond while loop
    ih = 3     # horizons are listed shallowest-to-deepest, first horizon is surface elevation
    while ih < nh   
        println("ih = $ih")
        # pull out this horizon only
        inds = findall(x -> typeof(x) == Float64, H[:,ih])
        h = H[inds,ih:ih+1]
        # x-location of the well (assume it's not deviated)
        x_well = mean(wellLog[:,1])
        # find the nearest two horizon nodes to the well 
        a = findmin(abs.(h[:,1] .- x_well))
        if x_well < h[a[2],1]
            ind1 = a[2] - 1; ind2 = a[2];
        elseif x_well > h[a[2],1]
            ind1 = a[2]; ind2 = a[2] + 1;
        elseif x_well == h[a[2],1]
            ind1 = a[2]; ind2 = a[2];
        end
        println("h: $h")
        # find depth to this horizon at location of well (assume it's not deviated for now)
        x1 = h[ind1,1]; x2 = h[ind2,1]; z1 = h[ind1,2]; z2 = h[ind2,2]; x = x_well;
        if ih == 3
            # shallowest horizon; zprev is the surface
            inds_ = findall(x -> typeof(x) == Float64, H[:,1])
            h_ = H[inds_,1:2]
            # find surface elevation at location of well (assume it's not deviated for now)
            z1_ = h_[ind1,2]; z2_ = h_[ind2,2];
            zprev = linearInterpolate(z1_,z2_,x1,x2,x)
            global zinterp = linearInterpolate(z1,z2,x1,x2,x)
        else
            zprev =  zinterp    # depth of previous horizon
            global zinterp = linearInterpolate(z1,z2,x1,x2,x)
        end
        # find all the well log values in this geologic unit
        p = findall(t -> t .< zinterp, wellLog[:,2])
        q = findall(t -> t .> zprev, wellLog[:,2])
        inds = intersect(p,q)
        println("Looking between $zprev and $zinterp, we have $(length(inds)) well log values")
        println("well log here: $(wellLog[:,2])")
        tmp = wellLog[inds,:]
        # tack an upper and lower row to represent the bounds of this geologic unit
        thisUnit = [tmp[1,1] zprev tmp[1,3:4]'; tmp; tmp[end,1] zinterp tmp[end,3:4]']
        # add this well log section to our vector of arrays
        push!(wellLogUnits,thisUnit)

        ih = ih + 2

    end

    # println("final geologic unit...")
    # # do it one last time for the last geologic unit (which has no bounding horizon on the bottom)
    # zprev = zinterp
    # zinterp = zmax
    # # find all the well log values in this geologic unit
    # p = findall(t -> t .< zinterp, wellLog[:,2])
    # q = findall(t -> t .> zprev, wellLog[:,2])
    # inds = intersect(p,q)
    # tmp = wellLog[inds,:]
    # # tack an upper and lower row to represent the bounds of this geologic unit
    # thisUnit = [tmp[1,1] zprev tmp[1,3:4]'; tmp; tmp[end,1] zinterp tmp[end,3:4]']
    # # add this well log section to our vector of arrays
    # push!(wellLogUnits,thisUnit)

    return wellLogUnits

end