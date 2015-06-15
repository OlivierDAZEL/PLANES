function Q=derive_polynom_2D_x(P)


nx=size(P,1);

Mx=(ones(nx,1)*(0:nx-1));
Q=Mx.*P;
Q(:,1)=[];