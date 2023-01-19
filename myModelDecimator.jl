
# This function takes as input a model defined as as N x 3 array
#   [ x (m),  y (m),  z (m) ]
# where each row is a point (x,y,z) in model space. The function
# outputs an N/(d^2) x 3 array that is the original model decimated by N
# in each of the x and y directions. 

function myModelDecimator(model_in,dec_x,dec_y)
    # this function assumes that the each row increments by dx until it resets to some smaller
    # x value and begins incrementing again, but this time at a y value that is dy larger 

    model_out = zeros(1,size(model_in,2))

    cntr = 0
    ind1 = 1
    ind2 = findnext(diff(model_in[:,1],dims=1) .< 0,ind1)
    while typeof(ind2) != Nothing
        if mod(cntr,dec_y) == 0
            # println("decimating from $ind1:$ind2 of $(size(model_in,1))")
            model_out = [model_out; model_in[ind1:dec_x:ind2,:]]
            ind1 = ind2 + 1
            ind2 = findnext(diff(model_in[:,1],dims=1) .< 0,ind1)
        else
            ind1 = ind2 + 1
            ind2 = findnext(diff(model_in[:,1],dims=1) .< 0,ind1)
        end
        cntr = cntr + 1
    end

    if mod(cntr,dec_y) == 0
        # println("decimating from $ind1:$ind2 of $(size(model_in,1)), last")
        model_out = [model_out; model_in[ind1:dec_x:end,:]]
    end

    # get rid of that first row of zeros
    model_out = model_out[2:end,:]

    return model_out

end

function myModelDecimator(model_in,nx,ny,dec_x,dec_y)
    # this function assumes the model grid is perfectly regular in x and y, with nx and ny
    # values in each x and y row of the grid and nx*ny values in model_in in total

    Nin = size(model_in,1)
    nxout = length(collect(1:dec_x:nx))
    nyout = length(collect(1:dec_y:ny))
    Nout = nxout*nyout
    model_tmp = zeros(nxout*ny,size(model_in,2))
    model_out = zeros(Nout,size(model_in,2))

    # first decimate along x 
    # println(" ")
    println("decimating along x...")
    # println(" ")
    cntr = 1
    for im=1:nx:Nin
        if mod(im,1e3) == 0
            println("$im of $Nin")
        end
        indend = im + nx - 1
        inds = collect(im:dec_x:indend)
        # println("inds = $inds")
        iout1 = (cntr-1)*nxout + 1
        iout2 = iout1 + nxout - 1
        # println("iout1:iout2 = $(collect(iout1:iout2))")
        global model_tmp[iout1:iout2,:] = model_in[inds,:]
        cntr = cntr + 1
    end

    # next decimate along y
    # println(" ")
    println("decimating along y...")
    # println(" ")
    cntr = 1
    for im=1:nxout:nxout*ny
        if mod(im,1e3) == 0
            println("$im of $(nxout*ny)")
        end
        ind1 = (cntr-1)*ny + 1
        ind2 = ind1 + nxout - 1
        # println("im:im+nxout-1 = $im:$(im+nxout-1)")
        # println("ind1:ind2 = $ind1:$ind2")
        model_out[im:im+nxout-1] = model_tmp[ind1:ind2]
        cntr = cntr + 1
    end

    return model_out

end