function f=integrate_polynom(P,lx)

% int_{x=0}^{x=lx} f(x) dx

nx=size(P,2);

Mx=(1:nx);
f=sum(P.*((lx.^Mx)./Mx));


