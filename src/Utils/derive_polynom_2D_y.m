function Q=derive_polynom_2D_y(P)


ny=size(P,2);

My=(ones(ny,1)*(0:ny-1))';
Q=My.*P;

Q(1,:)=[];