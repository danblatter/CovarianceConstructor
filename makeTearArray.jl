

function makeTearArray(C,tear)

n = size(C,1)
T = zeros(Int64,n)

if ~(tear[1] == Inf)        # tear[1] == Inf signals there is no tear function
    for i=1:n
        if tear[1,1] <= C[i,1] <= tear[1,end]
            xi = findfirst(tear[:,1] .> C[i,1])
            if ~isnothing(xi) && xi > 1         # otherwise this tear doesn't apply to this model parameter
                m = (tear[2,xi] - tear[2,xi-1])/(tear[1,xi] - tear[1,xi-1]) #slope
                b = tear[2,xi] - m*tear[1,xi]                               # y-intercept
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