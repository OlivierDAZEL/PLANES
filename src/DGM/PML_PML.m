% PML_PML.m
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
%

c_1=centre_element(e_1,nodes,elements);
[tau_x_e1,tau_y_e1]=parameter_PML(element_label(e_1));
ne1=sqrt((nx/tau_x_e1)^2+(ny/tau_y_e1)^2);


c_2=centre_element(e_2,nodes,elements);
[tau_x_e2,tau_y_e2]=parameter_PML(element_label(e_2));
ne2=sqrt((nx/tau_x_e2)^2+(ny/tau_y_e2)^2);


Re1e1=(air.Z*tau_x_e2*tau_y_e2*ne2-air.Z*tau_x_e1*tau_y_e1*ne1)/(air.Z*tau_x_e2*tau_y_e2*ne2+air.Z*tau_x_e1*tau_y_e1*ne1);
Re1e2=(2.00D+00)*air.Z*tau_x_e2*tau_y_e2*ne2^2/(air.Z*ne1*tau_x_e2*tau_y_e2*ne2+air.Z*tau_x_e1*tau_y_e1*ne1^2);
Re2e1=(2.00D+00)*air.Z*tau_x_e1*tau_y_e1*ne1^2/(air.Z*tau_x_e2*tau_y_e2*ne2^2+air.Z*ne2*tau_x_e1*tau_y_e1*ne1);
Re2e2=(air.Z*tau_x_e1*tau_y_e1*ne1-air.Z*tau_x_e2*tau_y_e2*ne2)/(air.Z*tau_x_e2*tau_y_e2*ne2+air.Z*tau_x_e1*tau_y_e1*ne1);
 
 
F_e1e1(1,1)=-air.Z*(nx^2)*(1+Re1e1)/(2*tau_x_e1^2*ne1);
F_e1e1(1,2)=-air.Z*(nx*ny)*(1+Re1e1)/(2*tau_x_e1*tau_y_e1*ne1);
F_e1e1(1,3)=-(nx)*(1+Re1e1)/(2*tau_x_e1);
F_e1e1(2,1)=-air.Z*(nx*ny)*(1+Re1e1)/(2*tau_x_e1*tau_y_e1*ne1);
F_e1e1(2,2)=-air.Z*(ny^2)*(1+Re1e1)/(2*tau_y_e1^2*ne1);
F_e1e1(2,3)=-(ny)*(1+Re1e1)/(2*tau_y_e1);
F_e1e1(3,1)=(nx)*(Re1e1-1)/(2*tau_x_e1);
F_e1e1(3,2)=(ny)*(Re1e1-1)/(2*tau_y_e1);
F_e1e1(3,3)=ne1*(Re1e1-1)/(2*air.Z);
 
 
F_e2e2(1,1)=-air.Z*(nx^2)*(1+Re2e2)/(2*tau_x_e2^2*ne2);
F_e2e2(1,2)=-air.Z*(nx*ny)*(1+Re2e2)/(2*tau_x_e2*tau_y_e2*ne2);
F_e2e2(1,3)=(nx)*(1+Re2e2)/(2*tau_x_e2);
F_e2e2(2,1)=-air.Z*(nx*ny)*(1+Re2e2)/(2*tau_x_e2*tau_y_e2*ne2);
F_e2e2(2,2)=-air.Z*(ny^2)*(1+Re2e2)/(2*tau_y_e2^2*ne2);
F_e2e2(2,3)=(ny)*(1+Re2e2)/(2*tau_y_e2);
F_e2e2(3,1)=-(nx)*(Re2e2-1)/(2*tau_x_e2);
F_e2e2(3,2)=-(ny)*(Re2e2-1)/(2*tau_y_e2);
F_e2e2(3,3)=ne2*(Re2e2-1)/(2*air.Z);
 
 
F_e2e1(1,1)=air.Z*(nx^2)/(2*tau_x_e2*tau_x_e1*ne1^2);
F_e2e1(1,2)=air.Z*(nx*ny)/(2*tau_x_e2*tau_y_e1*ne1^2);
F_e2e1(1,3)=(nx)*air.Z/(2*air.Z*tau_x_e2*ne1);
F_e2e1(2,1)=air.Z*(nx*ny)/(2*tau_x_e1*tau_y_e2*ne1^2);
F_e2e1(2,2)=air.Z*(ny^2)/(2*tau_y_e2*tau_y_e1*ne1^2);
F_e2e1(2,3)=(ny)*air.Z/(2*air.Z*tau_y_e2*ne1);
F_e2e1(3,1)=ne2*(nx)/(2*tau_x_e1*ne1^2);
F_e2e1(3,2)=ne2*(ny)/(2*tau_y_e1*ne1^2);
F_e2e1(3,3)=ne2*1/(2*air.Z*ne1);
 
F_e2e1=Re2e1*F_e2e1*ne2;
 
 
F_e1e2(1,1)=air.Z*(nx^2)/(2*tau_x_e1*tau_x_e2*ne2^2);
F_e1e2(1,2)=air.Z*(nx*ny)/(2*tau_x_e1*tau_y_e2*ne2^2);
F_e1e2(1,3)=-(nx)*air.Z/(2*air.Z*tau_x_e1*ne2);
F_e1e2(2,1)=air.Z*(nx*ny)/(2*tau_x_e2*tau_y_e1*ne2^2);
F_e1e2(2,2)=air.Z*(ny^2)/(2*tau_y_e1*tau_y_e2*ne2^2);
F_e1e2(2,3)=-(ny)*air.Z/(2*air.Z*tau_y_e1*ne2);
F_e1e2(3,1)=-ne1*(nx)/(2*tau_x_e2*ne2^2);
F_e1e2(3,2)=-ne1*(ny)/(2*tau_y_e2*ne2^2);
F_e1e2(3,3)=ne1*1/(2*air.Z*ne2);
 
F_e1e2=Re1e2*F_e1e2*ne1;


nx=cos(vec_theta);
ny=sin(vec_theta);

Phi=Phi_fluid_vector(nx,ny,air.Z,Shift_fluid);

II=int_edge_2vectorielle(1j*k_air*[nx*tau_x_e2;ny*tau_y_e2],-1j*k_air*[nx*tau_x_e2;ny*tau_y_e2],a,b,[c_2 c_2]);
MM=kron(II,F_e2e2);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_2);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_2);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi*tau_x_e2*tau_y_e2;

        
II=int_edge_2vectorielle(1j*k_air*[nx*tau_x_e2;ny*tau_y_e2],-1j*k_air*[nx*tau_x_e1;ny*tau_y_e1],a,b,[c_2 c_1]);
MM=kron(II,F_e2e1);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_2);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_1);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi*tau_x_e2*tau_y_e2;

II=int_edge_2vectorielle(1j*k_air*[nx*tau_x_e1;ny*tau_y_e1],-1j*k_air*[nx*tau_x_e2;ny*tau_y_e2],a,b,[c_1 c_2]);
MM=kron(II,F_e1e2);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_1);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_2);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi*tau_x_e1*tau_y_e1;

II=int_edge_2vectorielle(1j*k_air*[nx*tau_x_e1;ny*tau_y_e1],-1j*k_air*[nx*tau_x_e1;ny*tau_y_e1],a,b,[c_1 c_1]);
MM=kron(II,F_e1e1);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_1);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_1);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi*tau_x_e1*tau_y_e1;

