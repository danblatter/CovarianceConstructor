

using SuiteSparse

function getCovSquareRoot(B,dense_flag)

n = size(B,1)

if dense_flag == true
    A = cholesky(B,perm=1:n)     # Cholesky factorization of sparse covariance matrix B without permutation to reduce infill
    L = sparse(A.L)
else
    A = cholesky(B)     # Cholesky factorization making use of permutation to reduce infill
    Lp = sparse(A.L)
    P = sparse(1:n,A.p,ones(n))           # permutation matrix for permuted sparse Cholesky decomposition
    L = P'*Lp
end

return L

end