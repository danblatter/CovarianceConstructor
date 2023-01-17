# This code takes as input the 2D surface for the top of the seal and the bottom of the reservoir for
# the King Island model, as well as an angle θ measured positive counterclockwise from the x-direction,
# and produces as output a slice through that 2D surface along that angle. In other words, the
# output is an nx2 array where n is the number of (x,y) coordinates along the θ-profile where
# the top seal and bottom reservoir surfaces were interpolated.

using DelimitedFiles, LinearAlgebra

# read in the 2D slice (at angle θ) through the 3D model
# the columns are defined as follows:
# 1. x' (the new x-axis along the θ-direction)
# 2. z (m)
# 3. rho (ohm-m), linear (not log units)
# 4. x (original x-axis, θ=0)
# 5. y (original y-axis, θ=90)
println("loading in our 2D model slice")
myslice = readdlm("KI_60.dat", ',', Float64)
topSealSlice = zeros(size(myslice,1),2)
outputFilename = "KI60_topSeal_botRes.txt"

# read in the top of seal surface
# the columns of this array are defined as follows:
# 1. x (original x-axis, θ=0)
# 2. y (original y-axis, θ=90)
# 3. z (m)
# 4. repeat of column 3; omit
println("loading in the top seal surface")
topSeal = readdlm("top_shale_grd.dat")
topSeal = topSeal[:,1:3]
dx = abs(topSeal[1,1] - topSeal[2,1])
dy = dx

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

xrange = dx    # range in (m) from each model parameter within which to search
yrange = dy
# for (im,mvec) in enumerate(eachrow(myslice[:,4:5]))
for im=1:100
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
    topSealSlice[im,:] = [myslice[im,1] topSeal[indsxy[min_out[2]],3]]
end

writedlm(outputFilename,topSealSlice)

