# This code takes as input the 2D surface for the top of the seal and the bottom of the reservoir for
# the King Island model, as well as an angle θ measured positive counterclockwise from the x-direction,
# and produces as output a slice through that 2D surface along that angle. In other words, the
# output is an nx2 array where n is the number of (x,y) coordinates along the θ-profile where
# the top seal and bottom reservoir surfaces were interpolated.

using DelimitedFiles, LinearAlgebra, Statistics
include("myModelDecimator.jl")

# read in the 2D slice (at angle θ) through the 3D model
# the columns are defined as follows:
# 1. x' (the new x-axis along the θ-direction)
# 2. z (m)
# 3. rho (ohm-m), linear (not log units)
# 4. x (original x-axis, θ=0)
# 5. y (original y-axis, θ=90)
println("loading in our 2D model slice")
myslice = readdlm("/Users/Dan/Documents/Research/CarbonSequestration/modelViewing/slicesToView/KIslices/KI_60.dat", ',', Float64)
tmpind = findnext(diff(myslice[:,1]) .< 0, 1)
myslice = myslice[1:tmpind,:]
topSealbaseResSlice = zeros(tmpind,4)
topSealFilename = "Surfaces/KI60_topSeal.txt"
baseResFilename = "Surfaces/KI60_baseRes.txt"

# read in the top of seal surface
# the columns of this array are defined as follows:
# 1. x (original x-axis, θ=0)
# 2. y (original y-axis, θ=90)
# 3. z (m)
# 4. repeat of column 3; omit
println("loading in the top seal surface")
# read in the top seal surface file
topSeal = readdlm("top_shale_grd.dat")
# omit redundant 4th column
topSeal = topSeal[:,1:3]
# decimate to make the array manageable computationally
println("decimating the top of seal file")
topSeal = myModelDecimator(topSeal,5,5)

# read in the base of reservoir surface
# the columns of this array are defined as follows:
# 1. x (original x-axis, θ=0)
# 2. y (original y-axis, θ=90)
# 3. z (m)
# 4. repeat of column 3; omit
println("loading in the base reservoir surface")
# read in the base of reservoir surface file 
baseRes = readdlm("base_reservoir_grd.dat")
# remove redundant 4th column 
baseRes = baseRes[:,1:3]
# decimate the surface file to make it computationally manageable 
println("decimating the base of reservoir file")
baseRes = myModelDecimator(baseRes,5,5)

println("extracting top seal surface along 2D slice")

# # This version is so slow I don't even know if it works. Find a better way!
# for (im,mvec) in enumerate(eachrow(myslice[:,4:5]))
#     println("$im of $(size(myslice,1))")
#     if mod(im,1000) == 0
#         println("$im of $(size(myslice,1))")
#     end
#     for (js,svec) in enumerate(eachrow(topSeal[:,1:2]))
#         r[js] = norm(svec .- mvec)
#     end
#     tmp = findmin(abs.(r))
#     topSealSlice[im] = topSeal[tmp[2],3]
# end

xrange = median(diff(topSeal[:,1]))    # range in (m) from each model parameter within which to search
tmpind = findnext(diff(topSeal[:,2]) .> 0, 1)
yrange = abs(topSeal[tmpind+1,2] - topSeal[tmpind,2])
for (im,mvec) in enumerate(eachrow(myslice[:,4:5]))
# for im=1:100
    mvec = myslice[im,4:5]
    if mod(im,10) == 0
        println("$im of $(size(myslice,1))")
    end
    indnearx = findall(t -> abs.(mvec[1]-t) < xrange, topSeal[:,1])
    indneary = findall(t -> abs.(mvec[2]-t) < yrange, topSeal[:,2])
    indsxy = intersect(indnearx,indneary)
    tmpnodes = topSeal[indsxy,1:2]
    global r = zeros(length(indsxy))
    for (is,svec) in enumerate(eachrow(tmpnodes))
        global r[is] = norm(svec .- mvec)
    end
    min_out = findmin(abs.(r))
    # println("indsxy: $indsxy, min_out: $(indsxy[min_out[2]])")
    topSealbaseResSlice[im,1:2] = [myslice[im,1] topSeal[indsxy[min_out[2]],3]]
end

println("extracting base reservoir surface along 2D slice")

xrange = median(diff(baseRes[:,1]))    # range in (m) from each model parameter within which to search
tmpind = findnext(diff(baseRes[:,2]) .> 0, 1)
yrange = abs(baseRes[tmpind+1,2] - baseRes[tmpind,2])
for (im,mvec) in enumerate(eachrow(myslice[:,4:5]))
# for im=1:100
    mvec = myslice[im,4:5]
    if mod(im,10) == 0
        println("$im of $(size(myslice,1))")
    end
    indnearx = findall(t -> abs.(mvec[1]-t) < xrange, baseRes[:,1])
    indneary = findall(t -> abs.(mvec[2]-t) < yrange, baseRes[:,2])
    indsxy = intersect(indnearx,indneary)
    tmpnodes = baseRes[indsxy,1:2]
    global r = zeros(length(indsxy))
    for (is,svec) in enumerate(eachrow(tmpnodes))
        global r[is] = norm(svec .- mvec)
    end
    min_out = findmin(abs.(r))
    # println("indsxy: $indsxy, min_out: $(indsxy[min_out[2]])")
    topSealbaseResSlice[im,3:4] = [myslice[im,1] baseRes[indsxy[min_out[2]],3]]
end

# MARE2DEM assumes z is positive down, which is the opposite convention from the top seal and
# base reservoir surface files 
topSealbaseResSlice[:,2] = -topSealbaseResSlice[:,2]
topSealbaseResSlice[:,4] = -topSealbaseResSlice[:,4]

writedlm(topSealFilename,topSealbaseResSlice[:,1:2])
writedlm(baseResFilename,topSealbaseResSlice[:,3:4])

