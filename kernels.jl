# This file contains the various correlation kernels to use within buildCovariance

function GaspariCohn(r::Float64,l::Float64)
    # Gaspari-Cohn correlation kernel. Takes correlation length (l) and distance between two pts (r) as inputs
    if r < l
        b = -(1/4)*(r/l)^5 + (1/2)*(r/l)^4 +(5/8)*(r/l)^3 - (5/3)*(r/l)^2 + 1
    elseif r < 2*l
        b = (1/12)*(r/l)^5 - (1/2)*(r/l)^4 +(5/8)*(r/l)^3 + (5/3)*(r/l)^2 - 5*(r/l) + 4 - (2/3)*(l/r)
    else
        b = 0
    end
    return b
end

function GaspariCohn(R::Array{Float64,1},L::Array{Float64,1})
    # Gaspari-Cohn correlation kernel. Takes correlation length (l) and distance between two pts (r) as inputs
    if length(R) != length(L)
        println("if you want anisotropic covariance, the number of correlation lengths must equal")
        println("the number of spatial dimensions! Here the former is $(length(L)) while the latter is $(length(R))")
    end
    
    b = zeros(length(R))
    for ix in eachindex(R)
        r = R[ix]; l = L[ix]
        if r < l
            b[ix] = -(1/4)*(r/l)^5 + (1/2)*(r/l)^4 +(5/8)*(r/l)^3 - (5/3)*(r/l)^2 + 1
        elseif r < 2*l
            b[ix] = (1/12)*(r/l)^5 - (1/2)*(r/l)^4 +(5/8)*(r/l)^3 + (5/3)*(r/l)^2 - 5*(r/l) + 4 - (2/3)*(l/r)
        else
            b[ix] = 0
        end
    end

    return prod(b)

end

function MattiSpecial(r,li,lj,si,sj)
    # non-stationary covariance kernel; each model location has its own correlation length and variance
    prefactor = abs(li)^(0.25)*abs(lj)^(0.25)*abs((li+lj)/2)^(-0.5)
    expfactor = abs(r/((li+lj)/2))^(1.99)
    # expfactor = r' * r /((li+lj)/2)
    b = prefactor * exp(-expfactor) * si * sj
    return b
end

function squaredExponential(r,l)
    # squared exponential correlation kernel
    b = exp(-abs(r/l)^(1.99))
    return b
end

function exponential(r,l)
    # squared exponential correlation kernel
    b = exp(-abs((r/l)^(1.0)))
    return b
end

function buildGaspariCohn(C::Array{Float64,2},l::Array{Float64,1})

    println("lets compute using the Gaspari-Cohn kernel!")

    # check that the number of spatial dimensions is equal to the number of correlation lengths
    if size(C,2) != length(l)
        println("the number of spatial dimensions is not equal to the number of correlation lengths")
        println("# of spatial dimensions = $(size(C,2)); # of correlation lengths = $(length(l))")
    end

    n = size(C,1)
    M = zeros(Int64,n*n,1)
    N = zeros(Int64,n*n,1)
    V = zeros(n*n,1)

    k = 0
    for i=1:n
        if mod(i,100) == 0          # record our progress (this can take a while...)
            println("$i of $n")
        end
        for j=1:n
            # r = norm(C[i,:] - C[j,:])
            r = abs.(C[i,:] .- C[j,:])
            c = GaspariCohn(r,l)
            if c > 0                # only save non-zeros, since B is sparse
                k = k + 1
                M[k] = Int(i)
                N[k] = Int(j)
                V[k] = c
            end
        end
    end

    M = M[1:k]
    N = N[1:k]
    V = V[1:k]

    sparseratio = size(M,1)/(n^2)
    println("This sparse covariance takes $(100*sparseratio)% of the memory of the dense matrix")

    B = sparse(M,N,V,n,n)

    return B
end

function buildMattiSpecial(C,l::Number,s::Float64)
    println("lets compute using the Matti Special-GC kernel!")

    n = size(C,1)
    M = zeros(Int64,n*n,1)
    N = zeros(Int64,n*n,1)
    V = zeros(n*n,1)

    k = 0
    for i=1:n
        if mod(i,100) == 0          # record our progress (this can take a while...)
            println("$i of $n")
        end
        li = l          # correlation length for this model parameter location
        for j=1:n
            lj = l      # correlation length for this model parameter location
            r = norm(C[i,:] - C[j,:])   # distance between these two model parameters
            c = MattiSpecial(r,li,lj,s,s)
            cGC = GaspariCohn(r,l)
            c = c * cGC

            if c > 0                # only save non-zeros, since B is sparse
                k = k + 1
                M[k] = Int(i)
                N[k] = Int(j)
                V[k] = c
            else
                # println("model paramters $i and $j have exactly 0.0 covariance")
            end
        end
    end

    M = M[1:k]
    N = N[1:k]
    V = V[1:k]

    sparseratio = size(M,1)/(n^2)
    println("This sparse covariance takes $(100*sparseratio)% of the memory of the dense matrix")

    B = sparse(M,N,V,n,n)

    return B
end

function buildMattiSpecial(C,l::Function,s::Float64)
    println("lets compute using the Matti Special-GC kernel!")

    n = size(C,1)
    M = zeros(Int64,n*n,1)
    N = zeros(Int64,n*n,1)
    V = zeros(n*n,1)

    k = 0
    for i=1:n
        if mod(i,100) == 0          # record our progress (this can take a while...)
            println("$i of $n")
        end
        for j=1:n
            x = C[i,1]; z = C[i,2]; li = l(x,z)          # correlation length for this model parameter location
            x = C[j,1]; z = C[j,2]; lj = l(x,z)      # correlation length for this model parameter location
            r = norm(C[i,:] - C[j,:])   # distance between these two model parameters
            c = MattiSpecial(r,li,lj,s,s)
            lgc = (li+lj)/2; cGC = GaspariCohn(r,lgc)
            c = c * cGC
            
            if c > 0                # only save non-zeros, since B is sparse
                k = k + 1
                M[k] = Int(i)
                N[k] = Int(j)
                V[k] = c
            else
                # println("model paramters $i and $j have exactly 0.0 covariance")
            end
        end
    end

    M = M[1:k]
    N = N[1:k]
    V = V[1:k]

    sparseratio = size(M,1)/(n^2)
    println("This sparse covariance takes $(100*sparseratio)% of the memory of the dense matrix")

    B = sparse(M,N,V,n,n)

    return B
end

function buildMattiSpecial(C,l::Array{Float64,1},s::Array{Float64,1})
    println("lets compute using the Matti Special-GC kernel!")

    n = size(C,1)
    println("mean of std: $(mean(s))")
    M = zeros(Int64,n*n,1)
    N = zeros(Int64,n*n,1)
    V = zeros(n*n,1)

    k = 0
    for i=1:n
        # if mod(i,100) == 0          # record our progress (this can take a while...)
        #     println("$i of $n")
        # end
        for j=1:n
            li = l[i]          # correlation length for ith model parameter location
            lj = l[j]          # correlation length for jth model parameter location
            r = norm(C[i,:] - C[j,:])   # distance between these two model parameters
            si = s[i]
            sj = s[j]
            c = MattiSpecial(r,li,lj,si,sj)
            lgc = (li+lj)/2; cGC = GaspariCohn(r,2*lgc)
            c = c * cGC
            
            if c > 0                # only save non-zeros, since B is sparse
                k = k + 1
                M[k] = Int(i)
                N[k] = Int(j)
                V[k] = c
            else
                # println("model paramters $i and $j have exactly 0.0 covariance")
            end
        end
    end

    M = M[1:k]
    N = N[1:k]
    V = V[1:k]

    sparseratio = size(M,1)/(n^2)
    println("This sparse covariance takes $(100*sparseratio)% of the memory of the dense matrix")

    B = sparse(M,N,V,n,n)

    return B
end

function buildExponential(C,l)
    println("lets compute using the exponential kernel!")

    n = size(C,1)
    M = zeros(Int64,n*n,1)
    N = zeros(Int64,n*n,1)
    V = zeros(n*n,1)

    k = 0
    for i=1:n
        if mod(i,100) == 0          # record our progress (this can take a while...)
            println("$i of $n")
        end
        for j=1:n
            r = norm(C[i,:] - C[j,:])   # distance between these two model parameters
            c = exponential(r,l)
            if c > 0                    # only save non-zeros, since B is sparse
                k = k + 1
                M[k] = Int(i)
                N[k] = Int(j)
                V[k] = c
            else
                # println("model paramters $i and $j have exactly 0.0 covariance")
            end
        end
    end

    M = M[1:k]
    N = N[1:k]
    V = V[1:k]

    sparseratio = size(M,1)/(n^2)
    println("This sparse covariance takes $(100*sparseratio)% of the memory of the dense matrix")

    B = sparse(M,N,V,n,n)

    return B
end

function buildSquaredExponential(C,l)
    println("lets compute using the squared exponential kernel!")

    n = size(C,1)
    M = zeros(Int64,n*n,1)
    N = zeros(Int64,n*n,1)
    V = zeros(n*n,1)

    k = 0
    for i=1:n
        if mod(i,100) == 0          # record our progress (this can take a while...)
            println("$i of $n")
        end
        for j=1:n
            r = norm(C[i,:] - C[j,:])   # distance between these two model parameters
            c = squaredExponential(r,l)
            if c > 0                    # only save non-zeros, since B is sparse
                k = k + 1
                M[k] = Int(i)
                N[k] = Int(j)
                V[k] = c
            else
                # println("model paramters $i and $j have exactly 0.0 covariance")
            end
        end
    end

    M = M[1:k]
    N = N[1:k]
    V = V[1:k]

    sparseratio = size(M,1)/(n^2)
    println("This sparse covariance takes $(100*sparseratio)% of the memory of the dense matrix")

    B = sparse(M,N,V,n,n)

    return B
end