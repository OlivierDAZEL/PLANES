function f=integrate_polynom_2D(P,lx,ly)


nx=size(P,2);
ny=size(P,1);

Mx=(ones(ny,1)*(1:nx));
My=(ones(nx,1)*(1:ny))';


f=sum(sum(P.*((lx.^Mx)./Mx).*((ly.^My)./My)));
