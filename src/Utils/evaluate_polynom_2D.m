function f=evaluate_polynom_2D(P,x,y)


nx=size(P,2);
ny=size(P,1);

Mx=(ones(ny,1)*(0:nx-1));
My=(ones(nx,1)*(0:ny-1))';


f=sum(sum(P.*((x.^Mx).*(y.^My))));
