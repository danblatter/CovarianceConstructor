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
topSealSlice = zeros(size(myslice,1))

# read in the top of seal surface
# the columns of this array are defined as follows:
# 1. x (original x-axis, θ=0)
# 2. y (original y-axis, θ=90)
# 3. z (m)
# 4. repeat of column 3; omit
println("loading in the top seal surface")
topSeal = readdlm("top_shale_grd.dat")
topSeal = topSeal[:,1:3]

r = zeros(size(topSeal,1))

println("extracting top seal surface along 2D slice")
for (im,mvec) in enumerate(eachrow(myslice[:,4:5]))
    println("$im of $(size(myslice,1))")
    if mod(im,1000) == 0
        println("$im of $(size(myslice,1))")
    end
    for (js,svec) in enumerate(eachrow(topSeal[:,1:2]))
        r[js] = norm(svec .- mvec)
    end
    tmp = findmin(abs.(r))
    topSealSlice[im] = topSeal[tmp[2],3]
end

