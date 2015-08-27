function Q=derive_polynom(P)


nx=size(P,2);

Mx=(0:nx-1);
Q=Mx.*P;
Q(:,1)=[];