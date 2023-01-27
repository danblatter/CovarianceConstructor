
# # This script creates the following synthetic 1D profiles of:
# # 1. mean resistivity (ohm-m)
# # 2. standard deviation of resistivity (ohm-m)

# using DelimitedFiles, Statistics

# # output filename
# # outfilename = "KIrhoWellLog_fullyMeshed.txt"
# outfilename = "KIrhoWellLog_justReservoir.txt"
# # load the model parameter centroids (x,z)
# # C = readdlm("KI_60centroids_fullyMeshed.txt")
# C = readdlm("KI_60centroids_justReservoir.txt")
# # load two simple horizons (top and bottom of seal)
# # Rho = readdlm("KI_60rho_fullyMeshed.txt")
# Rho = readdlm("KI_60rho_justReservoir.txt")

# # load in the resistivity well log Erika gave me
# tmp = readdlm("ki128-res-log.dat")
# logρ_well = Float64.(tmp[2:end,2:3])
# logρ_well[:,1] = log10.(logρ_well[:,1])     # switch to log10 resistivity

# # define the depth range of this well log
# wellX = 0       # (m) we're not doing deviated wells right now...
# xRange = 125        # (m) make this at least as large as your grid spacing in the x-direction
# # row indices of all model parameter within xRange of the well position
# wellInds = findall(t -> abs.(t-wellX) < xRange, C[:,1])
# wellZs = C[wellInds,2]  # depth values of all model parameters within dx of the well position
# dz = 10         # (m) well log point density; make this less than the model grid spacing in z-direction
# z0 = minimum(wellZs) - dz/2
# zend = maximum(wellZs) + dz/2

# # (x,z) coordinates of well log values
# Z = Float64.(collect(z0:dz:zend))
# X = wellX .* ones(size(Z))
# nz = length(Z)
# # mean of log-resistivity
# meanRho = zeros(nz)
# # standard deviation of log-resistivity
# stdRho = zeros(nz)

# zRange = 50

# for iz=1:nz
#     if mod(iz,100) == 0
#         println("computing $iz of $nz")
#     end
#     local z = Z[iz]
#     # find all the model paramters within zRange of the well
#     local p = findall(t -> abs.(t-z) < zRange, C[:,2])
#     local inds = intersect(p,wellInds)
#     println("number used to compute this well log point: $(length(inds))")
#     meanRho[iz] = mean(Rho[inds])
# end

# # find the well log resistivity values that overlap with the synthetic well log
# inds1 = findall(t -> t > z0, logρ_well[:,2])
# inds2 = findall(t -> t < zend, logρ_well[:,2])
# indsoverlap = intersect(inds1,inds2)
# stdRho = stdRho .+ std(logρ_well[indsoverlap,1])

# rhoWellLog = [X Z meanRho stdRho]

# writedlm(outfilename,rhoWellLog)

function makeSyntheticKIwell(C,Rho,wellX)

    # In this function we are creating a 1D depth profile of resistivity (the mean) and uncertainty in resistivity 
    # (the standard deviation) over the depth range of the invertible model parameters at a particular lateral 
    # (surface) position. The purpose is to subsequently use this, in conjuction with geologic surfaces, to 
    # generate mean and standard deviation models defined at each point in the subsurface.
    #
    # In practice, one would simply use a resistivity well log to produce both the mean and standard deviation.
    # Here, we are using a vertical profile through the true model and a short segment of a resistivity log
    # (I don't know where the measurements were made, but I think it is synthetic in any case)

    # load in the resistivity well log Erika gave me
    tmp = readdlm("ki128-res-log.dat")
    logρ_well = Float64.(tmp[2:end,2:3])
    logρ_well[:,1] = log10.(logρ_well[:,1])     # switch to log10 resistivity

    xRange = 125        # (m) make this at least as large as your grid spacing in the x-direction
    zRange = 50         # (m) make this at least as large as your grid spacing in depth
    # row indices of all model parameter within xRange of the well position
    wellInds = findall(t -> abs.(t-wellX) < xRange, C[:,1])
    wellZs = C[wellInds,2]  # depth values of all model parameters within dx of the well position
    dz = 10         # (m) well log point density; make this less than the model grid spacing in z-direction
    z0 = minimum(wellZs) - dz/2     # make sure these include the entire model inversion region
    zend = maximum(wellZs) + dz/2

    # (x,z) coordinates of well log values
    Z = Float64.(collect(z0:dz:zend))
    X = wellX .* ones(size(Z))
    nz = length(Z)
    # mean of log-resistivity
    meanRho = zeros(nz)
    # standard deviation of log-resistivity
    stdRho = zeros(nz)

    # use the true model values nearest the well location (wellX) as a proxy for the prior model mean
    for iz=1:nz
        if mod(iz,100) == 0
            println("computing $iz of $nz")
        end
        local z = Z[iz]
        # find all the model paramters within zRange of the well
        local p = findall(t -> abs.(t-z) < zRange, C[:,2])
        local inds = intersect(p,wellInds)
        println("number used to compute this well log point: $(length(inds))")
        meanRho[iz] = mean(Rho[inds])
    end

    # use the sample standard deviation of the resistivity log over the depth range of the 
    # inversion model as a proxy for the prior model standard deviation
    # find the well log resistivity values that overlap with the synthetic well log
    inds1 = findall(t -> t > z0, logρ_well[:,2])
    inds2 = findall(t -> t < zend, logρ_well[:,2])
    indsoverlap = intersect(inds1,inds2)
    stdRho = stdRho .+ std(logρ_well[indsoverlap,1])

    rhoWellLog = [X Z meanRho stdRho]

    return rhoWellLog

end