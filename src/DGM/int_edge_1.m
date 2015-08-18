function f=int_edge_1(k_1,a,b,centre)


a_1=a-centre(1:2);

t=(b-a)/norm(b-a);

temp=(k_1).'*t;

if ((temp==0)|(norm(temp)<1e-6*(norm(k_1))))
    f=norm(b-a);
else
    f=(exp(norm(b-a)*temp)-1)/temp;
end

f=f*exp(k_1.'*a_1);

