% boundary_normal_velocity.m
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
%%%%% Coordinates of the edge's vertices

a=nodes(edges.dirichlets(ie,1),1:2)';
b=nodes(edges.dirichlets(ie,2),1:2)';

%%%%% Element linked to the edge

num_element=edges.dirichlets(ie,3);
center_element=mean(nodes(nonzeros(elem.nodes(num_element,:)),1:2))';

%%%%% Physical parameters

parameter_element

% (nx,ny) normal vector 

[nx,ny]=normal_edge_out_element(a,b,center_element);


P_e_in=[nx;ny;air.Z];
Q_e_in=[nx/2 ny/2 1/(2*air.Z)];
P_e_out=[nx;ny;-air.Z];
Q_e_out=[nx/2 ny/2 -1/(2*air.Z)];

C_e=[nx ny 0];

temp=inv(C_e*P_e_in);
R_e_tilde=-temp*C_e*P_e_in;
P_e_tilde=(P_e_in+P_e_out*R_e_tilde)*Q_e_in;

A_x=[0 0 1/rho_e;0 0 0; rho_e*c_e^2 0 0];
A_y=[0 0 0;0 0 1/rho_e; 0 rho_e*c_e^2 0];
F_e=(A_x*nx+A_y*ny);

F_e_tilde=M_e*F_e*P_e_tilde;


nx=cos(vec_theta);
ny=sin(vec_theta);
Phi=Phi_fluid_vector(nx,ny,Z_e,Shift_fluid);

II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a,b,[center_element center_element]);
MM=kron(II,F_e_tilde);
indice_test  =((1:data_model.theta_DGM.nb)-1)+dof_start_element(num_element);  
indice_champs=((1:data_model.theta_DGM.nb)-1)+dof_start_element(num_element);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;

