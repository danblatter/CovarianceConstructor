
function linearInterpolate(y0,y1,x0,x1,x)

    if x0 == x1
        y = y0
    else
        y = y0 + (x-x0)*(y1-y0)/(x1-x0)
    end

    return y

end