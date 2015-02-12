function [F_plus,F_moins]=Split_PML(nx,ny,Z_e,tau_x,tau_y)

ne=sqrt((nx/tau_x)^2+(ny/tau_y)^2);


F_plus(1,1)=-Z_e*(nx^2/2)/(ne*tau_x^2);
F_plus(1,2)=-Z_e*(nx*ny/2)/(ne*tau_x*tau_y);
F_plus(1,3)=(nx/2)/tau_x;
F_plus(2,1)=-Z_e*(nx*ny/2)/(ne*tau_x*tau_y);
F_plus(2,2)=-Z_e*(ny^2/2)/(ne*tau_y^2);
F_plus(2,3)=(ny/2)/tau_y;
F_plus(3,1)=(nx/2)/tau_x;
F_plus(3,2)=(ny/2)/tau_y;
F_plus(3,3)=-ne*(0.5)/Z_e;


F_moins(1,1)=Z_e*(nx^2/2)/(ne*tau_x^2);
F_moins(1,2)=Z_e*(nx*ny/2)/(ne*tau_x*tau_y);
F_moins(1,3)=(nx/2)/tau_x;
F_moins(2,1)=Z_e*(nx*ny/2)/(ne*tau_x*tau_y);
F_moins(2,2)=Z_e*(ny^2/2)/(ne*tau_y^2);
F_moins(2,3)=(ny/2)/tau_y;
F_moins(3,1)=(nx/2)/tau_x;
F_moins(3,2)=(ny/2)/tau_y;
F_moins(3,3)=ne*(0.5)/Z_e;






