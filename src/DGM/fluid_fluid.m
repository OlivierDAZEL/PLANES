% fluid_fluid.m
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

e_edge=e_1;
c_1=centre_element(e_1,nodes,elements);

parameter_element

Z_e1=Z_e;
k_e1=k_e;

e_edge=e_2;
c_2=centre_element(e_2,nodes,elements);
parameter_element

Z_e2=Z_e;
k_e2=k_e;


unc=1;
deuxc=2;


Re1e1=(Z_e2-Z_e1)/(Z_e1+Z_e2);
Re1e2=(2)*Z_e2/(Z_e1+Z_e2);
Re2e1=(2)*Z_e1/(Z_e1+Z_e2);
Re2e2=(Z_e1-Z_e2)/(Z_e1+Z_e2);



F_e1e1(1,1)=-Z_e1*(nx^2)*(unc+Re1e1)/deuxc;
F_e1e1(1,2)=-Z_e1*(nx*ny)*(unc+Re1e1)/deuxc;
F_e1e1(1,3)=-(nx)*(unc+Re1e1)/deuxc;
F_e1e1(2,1)=-Z_e1*(nx*ny)*(unc+Re1e1)/deuxc;
F_e1e1(2,2)=-Z_e1*(ny^2)*(unc+Re1e1)/deuxc;
F_e1e1(2,3)=-(ny)*(unc+Re1e1)/deuxc;
F_e1e1(3,1)=(nx)*(Re1e1-unc)/deuxc;
F_e1e1(3,2)=(ny)*(Re1e1-unc)/deuxc;
F_e1e1(3,3)=(Re1e1-unc)/(deuxc*Z_e1);


F_e2e2(1,1)=-Z_e2*(nx^2)*(unc+Re2e2)/deuxc;
F_e2e2(1,2)=-Z_e2*(nx*ny)*(unc+Re2e2)/deuxc;
F_e2e2(1,3)=(nx)*(unc+Re2e2)/deuxc;
F_e2e2(2,1)=-Z_e2*(nx*ny)*(unc+Re2e2)/deuxc;
F_e2e2(2,2)=-Z_e2*(ny^2)*(unc+Re2e2)/deuxc;
F_e2e2(2,3)=(ny)*(unc+Re2e2)/deuxc;
F_e2e2(3,1)=-(nx)*(Re2e2-unc)/deuxc;
F_e2e2(3,2)=-(ny)*(Re2e2-unc)/deuxc;
F_e2e2(3,3)=(Re2e2-unc)/(deuxc*Z_e2);


F_e2e1(1,1)=Z_e2*(nx^2)/deuxc;
F_e2e1(1,2)=Z_e2*(nx*ny)/deuxc;
F_e2e1(1,3)=(nx)*Z_e2/(deuxc*Z_e1);
F_e2e1(2,1)=Z_e2*(nx*ny)/deuxc;
F_e2e1(2,2)=Z_e2*(ny^2)/deuxc;
F_e2e1(2,3)=(ny)*Z_e2/(deuxc*Z_e1);
F_e2e1(3,1)=(nx)/deuxc;
F_e2e1(3,2)=(ny)/deuxc;
F_e2e1(3,3)=unc/(deuxc*Z_e1);

F_e2e1=Re2e1*F_e2e1;


F_e1e2(1,1)=Z_e1*(nx^2)/deuxc;
F_e1e2(1,2)=Z_e1*(nx*ny)/deuxc;
F_e1e2(1,3)=-(nx)*Z_e1/(deuxc*Z_e2);
F_e1e2(2,1)=Z_e1*(nx*ny)/deuxc;
F_e1e2(2,2)=Z_e1*(ny^2)/deuxc;
F_e1e2(2,3)=-(ny)*Z_e1/(deuxc*Z_e2);
F_e1e2(3,1)=-(nx)/deuxc;
F_e1e2(3,2)=-(ny)/deuxc;
F_e1e2(3,3)=unc/(deuxc*Z_e2);

F_e1e2=Re1e2*F_e1e2;

nx=cos(vec_theta);
ny=sin(vec_theta);

Phi_1=Phi_fluid_vector(nx,ny,Z_e1,Shift_fluid);
Phi_2=Phi_fluid_vector(nx,ny,Z_e2,Shift_fluid);

II=int_edge_2vectorielle(1j*k_e2*[nx;ny],-1j*k_e2*[nx;ny],a,b,[c_2 c_2]);
MM=kron(II,F_e2e2);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_2);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_2);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi_2.'*MM*Phi_2;

        
II=int_edge_2vectorielle(1j*k_e2*[nx;ny],-1j*k_e1*[nx;ny],a,b,[c_2 c_1]);
MM=kron(II,F_e2e1);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_2);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_1);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi_2.'*MM*Phi_1;

II=int_edge_2vectorielle(1j*k_e1*[nx;ny],-1j*k_e2*[nx;ny],a,b,[c_1 c_2]);
MM=kron(II,F_e1e2);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_1);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_2);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi_1.'*MM*Phi_2;

II=int_edge_2vectorielle(1j*k_e1*[nx;ny],-1j*k_e1*[nx;ny],a,b,[c_1 c_1]);
MM=kron(II,F_e1e1);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_1);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_1);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi_1.'*MM*Phi_1;



