function [F_plus,F_moins]=Split_fluid(nx,ny,k_0,Z_0,omega,M)



k_splitting=diag([0 -omega/k_0 omega/k_0]);


Phi_splitting(1:3,1)=Phi_fluid_0(nx,ny);
Phi_splitting(1:3,2)=Phi_fluid( nx, ny,Z_0);
Phi_splitting(1:3,3)=Phi_fluid(-nx,-ny,Z_0);

Psi_splitting=inv(Phi_splitting);

F_plus =M*(-omega/k_0*Phi_splitting(:,2)*Psi_splitting(2,:));
F_moins=M*( omega/k_0*Phi_splitting(:,3)*Psi_splitting(3,:));







