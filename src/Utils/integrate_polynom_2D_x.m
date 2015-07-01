function f=integrate_polynom_2D_x(P,lx,y)


% int_{x=0}^{x=lx} f(x,y) dx


nx=size(P,2);
ny=size(P,1);

My=(ones(nx,1)*(0:ny-1))';


f=sum(P.*(y.^My));

Mx=(1:nx);
f=sum(f.*((lx.^Mx)./Mx));


