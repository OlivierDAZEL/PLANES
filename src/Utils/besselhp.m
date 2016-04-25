function f=besselhp(m,k,r)

f=(besselh(m-1,k,r)-besselh(m+1,k,r))/2;