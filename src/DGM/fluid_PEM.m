% fluid_PEM.m
%
% Copyright (C) 2015 < Olivier DAZEL <olivier.dazel@univ-lemans.fr> >
%
% This file is part of PLANES.
%
% PLANES (Porous LAum NumErical Simulator) is a software to compute the
% vibroacoustic response of sound packages containing coupled
% acoustic/elastic/porous substructures. It is mainly based on the
% Finite-Element Method and some numerical methods developped at
% LAUM (http://laum.univ-lemans.fr).
%
% You can download the latest version of PLANES at
%                 https://github.com/OlivierDAZEL/PLANES
% or find more details on Olivier's webpage
%                  http://perso.univ-lemans.fr/~odazel/
%
% For any questions or if you want to
% contribute to this project, contact Olivier.
%
% PLANES is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% e'=e_1
% e=e_2


e_edge=e_1;
parameter_element

M_e_p=M_e;

W_e_p_plus= Phi_fluid( nx, ny,Z_e);
W_e_p_moins=Phi_fluid(-nx,-ny,Z_e);
W_e_p_0=Phi_fluid_0(nx,ny);


Omega_e_p=inv([W_e_p_plus W_e_p_moins W_e_p_0]);
Omega_e_p_plus=Omega_e_p(1,:);
Omega_e_p_moins=Omega_e_p(2,:);

Lambda_e_p_moins=omega/k_e;
Lambda_e_p_plus=-Lambda_e_p_moins;

k_e1=k_e;
Z_e1=Z_e;

e_edge=e_2;
parameter_element


W_e_plus= Phi_Biot(-nx,-ny,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
W_e_moins=Phi_Biot( nx, ny,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
W_e_0=Phi_Biot_0(nx,ny);

Omega_e=inv([W_e_plus W_e_moins W_e_0]);
Omega_e_plus=Omega_e(1:3,:);
Omega_e_moins=Omega_e(4:6,:);

Lambda_e_moins=-diag([omega/delta_1 omega/delta_2 omega/delta_3]);
Lambda_e_plus=-Lambda_e_moins;

k_e=[delta_1 delta_2 delta_3];


B_e_p=zeros(4,3);
B_e=zeros(4,8);

B_e_p(3,1)=nx;
B_e_p(3,2)=ny;
B_e_p(4,3)=1;

B_e(1,1)=1;
B_e(2,2)=1;
B_e(3,3)=nx;
B_e(3,4)=ny;
B_e(4,8)=1;


M1=[-B_e*M_e*W_e_plus*Lambda_e_plus B_e_p*M_e_p*W_e_p_plus*Lambda_e_p_plus];
M2=[ B_e*M_e*W_e_moins*Lambda_e_moins -B_e_p*M_e_p*W_e_p_moins*Lambda_e_p_moins];

Refl=inv(M2)*M1;

R_ee=Refl(1:3,1:3);
R_eep=Refl(1:3,4);
R_epe=Refl(4,1:3);
R_epep=Refl(4,4);


F_e2e2=  M_e*(W_e_plus*Lambda_e_plus+W_e_moins*Lambda_e_moins*R_ee)*Omega_e_plus;
F_e2e1= M_e*W_e_moins*Lambda_e_moins*R_eep*Omega_e_p_plus;
F_e1e2=-M_e_p*W_e_p_moins*Lambda_e_p_moins*R_epe*Omega_e_plus;
F_e1e1=-M_e_p*(W_e_p_moins*Lambda_e_p_moins*R_epep+W_e_p_plus*Lambda_e_p_plus)*Omega_e_p_plus;

nx=cos(vec_theta);
ny=sin(vec_theta);

Phi_1=Phi_fluid_vector(nx,ny,Z_e1,Shift_fluid);
Phi_2=Phi_Biot_vector(nx,ny,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega,Shift_Biot);

II=int_edge_2vectorielle(1j*[delta_1*[nx;ny] delta_2*[nx;ny] delta_3*[nx;ny]],-1j*[delta_1*[nx;ny] delta_2*[nx;ny] delta_3*[nx;ny]],a,b,[c_2 c_2]);
MM=kron(II,F_e2e2);
indice_test  =[1+3*(0:nb_theta-1) 2+3*(0:nb_theta-1) 3+3*(0:nb_theta-1)]+dof_start_element(e_2)-1;
indice_champs=[1+3*(0:nb_theta-1) 2+3*(0:nb_theta-1) 3+3*(0:nb_theta-1)]+dof_start_element(e_2)-1;
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi_2.'*MM*Phi_2;

        
II=int_edge_2vectorielle(1j*[delta_1*[nx;ny] delta_2*[nx;ny] delta_3*[nx;ny]],-1j*k_e1*[nx;ny],a,b,[c_2 c_1]);
MM=kron(II,F_e2e1);
indice_test  =[1+3*(0:nb_theta-1) 2+3*(0:nb_theta-1) 3+3*(0:nb_theta-1)]+dof_start_element(e_2)-1; 
indice_champs=((1:nb_theta)-1)+dof_start_element(e_1);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi_2.'*MM*Phi_1;

II=int_edge_2vectorielle(1j*k_e1*[nx;ny],-1j*[delta_1*[nx;ny] delta_2*[nx;ny] delta_3*[nx;ny]],a,b,[c_1 c_2]);
MM=kron(II,F_e1e2);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_1);  
indice_champs=[1+3*(0:nb_theta-1) 2+3*(0:nb_theta-1) 3+3*(0:nb_theta-1)]+dof_start_element(e_2)-1;
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi_1.'*MM*Phi_2;

II=int_edge_2vectorielle(1j*k_e1*[nx;ny],-1j*k_e1*[nx;ny],a,b,[c_1 c_1]);
MM=kron(II,F_e1e1);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_1);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_1);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi_1.'*MM*Phi_1;
