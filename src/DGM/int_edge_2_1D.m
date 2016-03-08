function f=int_edge_2_1D(k_1,k_2,a,b,centres)

a_1=a-centres(1);
a_2=a-centres(2);

temp=(k_1+k_2);

if (norm(temp)<1e-6*(norm(k_1)))
    f=b-a;
else
    f=(exp(norm(b-a)*temp)-1)/temp;
end

f=f*exp(k_1.'*a_1)*exp(k_2.'*a_2);

