

function computeTearDistance(y::Vector{Float64},tear)
    d = Inf
    if tear[1,1] <= y[1] <= tear[end,1]
        xi = findfirst(tear[:,1] .> y[1])
        if ~isnothing(xi) && xi > 1         # otherwise this tear doesn't apply to this model parameter
            m = (tear[xi,2] - tear[xi-1,2])/(tear[xi,1] - tear[xi-1,1]) #slope
            b = tear[xi,2] - m*tear[xi,1]                               # y-intercept
            zinterp = m*y[1] + b          # tear surface depth at this model parameter
            d = abs(zinterp - y[2])
        end
    end
    return d
end


function makeTearArray(C,tear)

n = size(C,1)
T = zeros(Int64,n)

if ~(tear[1] == Inf)        # tear[1] == Inf signals there is no tear function
    for i=1:n
        if tear[1,1] <= C[i,1] <= tear[end,1]
            xi = findfirst(tear[:,1] .> C[i,1])
            if ~isnothing(xi) && xi > 1         # otherwise this tear doesn't apply to this model parameter
                m = (tear[xi,2] - tear[xi-1,2])/(tear[xi,1] - tear[xi-1,1]) #slope
                b = tear[xi,2] - m*tear[xi,1]                               # y-intercept
                zinterp = m*C[i,1] + b          # tear surface depth at this model parameter
                if C[i,2] >= zinterp
                    T[i] = 1
                else
                    T[i] = -1
                end
            end
        end
    end
end

return T

end