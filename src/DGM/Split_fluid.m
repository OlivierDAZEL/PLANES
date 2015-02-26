function [F_plus,F_moins]=Split_fluid(nx,ny,Z_e)


F_moins(1,1)=Z_e*(nx^2/2);
F_moins(1,2)=Z_e*(nx*ny/2);
F_moins(1,3)=(nx/2);
F_moins(2,1)=Z_e*(nx*ny/2);
F_moins(2,2)=Z_e*(ny^2/2);
F_moins(2,3)=(ny/2);
F_moins(3,1)=(nx/2);
F_moins(3,2)=(ny/2);
F_moins(3,3)=1/(2*Z_e);

F_plus(1,1)=-Z_e*(nx^2/2);
F_plus(1,2)=-Z_e*(nx*ny/2);
F_plus(1,3)= (nx/2);
F_plus(2,1)=-Z_e*(nx*ny/2);
F_plus(2,2)=-Z_e*(ny^2/2);
F_plus(2,3)= (ny/2);
F_plus(3,1)= (nx/2);
F_plus(3,2)= (ny/2);
F_plus(3,3)=-1/(2*Z_e);




