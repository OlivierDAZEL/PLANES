function Q=derive_polynom_2D_y_2(P)


ny=size(P,2);

My=(ones(ny,1)*(0:ny-1))';
Q=My.*P;
Q(end+1,:)=0;
Q(1,:)=[];