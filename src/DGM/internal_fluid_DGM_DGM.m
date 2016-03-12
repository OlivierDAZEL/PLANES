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
c_e=air.c;
k_e=omega/c_e;
Z_e=air.Z;
rho_e=air.rho;
M_e=diag([air.rho,air.rho,1/air.K]);



a=nodes(edges.flux(ie,1),1:2)';
b=nodes(edges.flux(ie,2),1:2)';


e_1=edges.flux(ie,3);
e_2=edges.flux(ie,4);
c_1=mean(nodes(nonzeros(elem.nodes(e_1,:)),1:2))';
c_2=mean(nodes(nonzeros(elem.nodes(e_2,:)),1:2))';
[nx,ny]=normal_edge_out_element(a,b,c_1);

M_1=M_e;
M_2=M_e;
C_2=[nx ny 0;0 0 1];
C_1=[nx ny 0;0 0 1];

P_1_in=[nx;ny;air.Z];
Q_1_in=[nx/2 ny/2 1/(2*air.Z)];
P_1_out=[nx;ny;-air.Z];
Q_1_out=[nx/2 ny/2 -1/(2*air.Z)];

P_2_in=P_1_out;
Q_2_in=Q_1_out;
P_2_out=P_1_in;
Q_2_out=Q_1_in;


A_x=[0 0 1/air.rho;0 0 0; air.K 0 0];
A_y=[0 0 0;0 0 1/air.rho; 0 air.K 0];

F_1= (A_x*nx+A_y*ny);
F_2=-(A_x*nx+A_y*ny);


P_11_tilde=P_1_in*Q_1_in;
P_12_tilde=P_1_out*Q_2_in;
P_21_tilde=P_2_out*Q_1_in;
P_22_tilde=P_2_in*Q_2_in;


F_11_tilde= M_1*F_1*P_11_tilde;
F_12_tilde= M_1*F_1*P_12_tilde;
F_21_tilde= M_2*F_2*P_21_tilde;
F_22_tilde= M_2*F_2*P_22_tilde;


nx=cos(vec_theta);
ny=sin(vec_theta);
Phi=Phi_fluid_vector(nx,ny,Z_e,Shift_fluid);

II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a,b,[c_2 c_2]);
MM=kron(II,F_22_tilde);
indice_test  =((1:data_model.theta_DGM.nb)-1)+dof_start_element(e_2);  
indice_champs=((1:data_model.theta_DGM.nb)-1)+dof_start_element(e_2);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;
        
II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a,b,[c_2 c_1]);
MM=kron(II,F_21_tilde);
indice_test  =((1:data_model.theta_DGM.nb)-1)+dof_start_element(e_2);  
indice_champs=((1:data_model.theta_DGM.nb)-1)+dof_start_element(e_1);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;

II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a,b,[c_1 c_2]);
MM=kron(II,F_12_tilde);
indice_test  =((1:data_model.theta_DGM.nb)-1)+dof_start_element(e_1);  
indice_champs=((1:data_model.theta_DGM.nb)-1)+dof_start_element(e_2);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;

II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a,b,[c_1 c_1]);
MM=kron(II,F_11_tilde);
indice_test  =((1:data_model.theta_DGM.nb)-1)+dof_start_element(e_1);  
indice_champs=((1:data_model.theta_DGM.nb)-1)+dof_start_element(e_1);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;
