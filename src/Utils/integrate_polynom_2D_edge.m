function f = integrate_polynom_2D_edge(P,coord_edge)

a=coord_edge(:,1);
b=coord_edge(:,2);
length_edge=norm(b-a);

nb_elem=100;
h=length_edge/nb_elem;
delta=(b-a)/nb_elem;

I_t=0;
I_rm=0;

for ii=1:nb_elem
    point=a+delta*(ii-1);
    I_t=I_t+evaluate_polynom_2D(P,point(1),point(2));
    point=a+delta*(ii  );
    I_t=I_t+evaluate_polynom_2D(P,point(1),point(2));
    point=a+delta*(ii-1/2);
    I_rm=I_rm+evaluate_polynom_2D(P,point(1),point(2));
end
I_t=I_t*h/2;
I_rm=I_rm*h;

f=(I_t+2*I_rm)/3;
end

