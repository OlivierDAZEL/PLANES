% internal_fluid.m
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
e_edge=e_1;
parameter_element

[F_moins,F_plus]=Split_fluid(nx,ny,Z_e);


C_e2=[0 1 0;0 0 1];
C_e1=[0 1 0;0 0 1];
W_e2=[nx -ny -nx;ny nx -ny;air.Z 0 air.Z];
Omega_e2=[nx/2 ny/2 1/(2*air.Z);-ny nx 0;-nx/2 -ny/2 1/(2*air.Z)];
W_e2_plus=W_e2(:,3);
W_e2_moins=W_e2(:,1);
Omega_e2_plus=Omega_e2(3,:);
Omega_e2_moins=Omega_e2(1,:);

W_e1=W_e2;
Omega_e1=Omega_e2;

W_e1_moins=W_e1(:,1);
W_e1_plus=W_e1(:,3);
Omega_e1_moins=Omega_e1(1,:);
Omega_e1_plus=Omega_e1(3,:);


% F_e=ny*[0 0 0 ;0 0 1/air.rho;0 air.K 0];
% 
% 
% dfgdgfdfggdf

nx=cos(vec_theta);
ny=sin(vec_theta);
Phi=Phi_fluid_vector(nx,ny,Z_e,Shift_fluid);

II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a,b,[c_2 c_2]);
MM=kron(II,F_plus);
indice_test  =((1:theta_DGM.nb)-1)+dof_start_element(e_2);  
indice_champs=((1:theta_DGM.nb)-1)+dof_start_element(e_2);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;
        
II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a,b,[c_2 c_1]);
MM=kron(II,F_moins);
indice_test  =((1:theta_DGM.nb)-1)+dof_start_element(e_2);  
indice_champs=((1:theta_DGM.nb)-1)+dof_start_element(e_1);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;

II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a,b,[c_1 c_2]);
MM=kron(II,-F_plus);
indice_test  =((1:theta_DGM.nb)-1)+dof_start_element(e_1);  
indice_champs=((1:theta_DGM.nb)-1)+dof_start_element(e_2);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;

II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a,b,[c_1 c_1]);
MM=kron(II,-F_moins);
indice_test  =((1:theta_DGM.nb)-1)+dof_start_element(e_1);  
indice_champs=((1:theta_DGM.nb)-1)+dof_start_element(e_1);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;
