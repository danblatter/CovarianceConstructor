# This file contains the various correlation kernels to use within buildCovariance

function GaspariCohn(r,l)
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

function MattiSpecial(r,li,lj)
    # non-stationary covariance kernel; each model location has its own correlation length and variance
    prefactor = abs(li)^(0.25)*abs(lj)^(0.25)*abs((li+lj)/2)^(-0.5)
    expfactor = abs(r/((li+lj)/2))^(1.99)
    # expfactor = r' * r /((li+lj)/2)
    b = prefactor * exp(-expfactor)
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

function buildGaspariCohn(C,l)

    println("lets compute using the Gaspari-Cohn kernel!")

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
            r = norm(C[i,:] - C[j,:])
            c = GaspariCohn(r,l)
#=             if abs(T[i] - T[j]) == 2    # these two model parameters are on opposite sides of a tear
                c = 0.0*c               # skip to next model parameter (correlation is zero by default)
            end =#
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

function buildMattiSpecial(C,l::Number,T)
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
            if abs(T[i] - T[j]) == 2                 # these two model parameters are on opposite sides of a tear
                li = 100; lj = 100
            end
            r = norm(C[i,:] - C[j,:])   # distance between these two model parameters
            c = MattiSpecial(r,li,lj)
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

function buildMattiSpecial(C,l::Function,T)
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
            if abs(T[i] - T[j]) == 2                 # these two model parameters are on opposite sides of a tear
                lj = 100; li = 100
            end
            r = norm(C[i,:] - C[j,:])   # distance between these two model parameters
            c = MattiSpecial(r,li,lj)
            lgc = (li+lj)/2; cGC = GaspariCohn(r,lgc)
            c = c * cGC

            if i == 3767 && j == 3635
                println("i=$i; j=$j; li=$(li); lj=$(lj); cGC=$(cGC); lgc=$(lgc); tear=$(abs(T[i] - T[j]) == 2)")
            end
            
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