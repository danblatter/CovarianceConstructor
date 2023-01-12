
# This script creates the following synthetic 1D profiles of:
# 1. mean resistivity (ohm-m)
# 2. standard deviation of resistivity (ohm-m)

using DelimitedFiles

# output filename
outfilename = "KIrhoWellLog.txt"
# load the model parameter centroids (x,z)
C = readdlm("KI_60centroids_fullyMeshed.txt")
# load two simple horizons (top and bottom of seal)
Rho = readdlm("KI_60rho_fullyMeshed.txt")

# (x,z) coordinates of well log values
Z = Float64.(collect(10:20:2200))
nz = length(Z)
X = zeros(nz)
# mean of log-resistivity
meanRho = zeros(nz)
# standard deviation of log-resistivity
stdRho = zeros(nz)

dx = 50
dz = 50

for iz=1:nz
    if mod(iz,100) == 0
        println("computing $iz of $nz")
    end
    local x = X[iz]
    local z = Z[iz]
    # find all the model paramters within dx of the well
    local q = findall(t -> abs.(t-x) < dx, C[:,1])
    local p = findall(t -> abs.(t-z) < dz, C[:,2])
    local inds = intersect(p,q)
    println("number used to compute this well log point: $(length(inds))")
    meanRho[iz] = mean(Rho[inds])
    stdRho[iz] = 0.1*(abs.(maximum(Rho) - minimum(Rho)))
end

rhoWellLog = [X Z meanRho stdRho]

writedlm(outfilename,rhoWellLog)