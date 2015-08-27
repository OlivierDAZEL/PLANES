function f=int_edge_2(k_1,k_2,a,b,centres)

a_1=a-centres(1:2,1);
a_2=a-centres(1:2,2);

t=(b-a)/norm(b-a);

temp=(k_1+k_2).'*t;

if (norm(temp)<1e-6*(norm(k_1)))
    f=norm(b-a);
else
    f=(exp(norm(b-a)*temp)-1)/temp;
end

f=f*exp(k_1.'*a_1)*exp(k_2.'*a_2);

