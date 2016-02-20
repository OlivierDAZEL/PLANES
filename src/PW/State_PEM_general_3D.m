% State_PEM_general_3D.m
%
% Copyright (C) 2014 < Olivier DAZEL <olivier.dazel@univ-lemans.fr> >
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
% https://github.com/OlivierDAZEL/PLANES
% or find more details on Olivier's webpage
% http://perso.univ-lemans.fr/~odazel/
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
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
% 
%
% The State vector for a fluid according to Dazel et al. JAP 2013
% S={u_z,p}
%
%
%


%%


function [k_z,SV]=State_PEM_general_3D(k_x,k_y,num_mat,omega,air)

eval(['Mat_porous_' num2str(num_mat)]);
properties_eqf
properties_PEM_aniso


M=  zeros(18,18);
A_x=zeros(18,18);
A_y=zeros(18,18);
A_z=zeros(18,18);

M(1:3,1:3)=-omega^2*rho_s_til;
M(1:3,4:6)=-omega^2*gamma_til*rho_eq_til;
M(4:6,1:3)=-omega^2*gamma_til*rho_eq_til;
M(4:6,4:6)=-omega^2*rho_eq_til;


A_x(1,7)=-1;
A_y(1,12)=-1;
A_z(1,11)=-1;

A_x(2,12)=-1;
A_y(2,8)=-1;
A_z(2,10)=-1;

A_x(3,11)=-1;
A_y(3,10)=-1;
A_z(3,9)=-1;

A_x(4,7+6)=-1;
A_y(4,12+6)=-1;
A_z(4,11+6)=-1;

A_x(5,12+6)=-1;
A_y(5,8+6)=-1;
A_z(5,10+6)=-1;

A_x(6,11+6)=-1;
A_y(6,10+6)=-1;
A_z(6,9+6)=-1;

M(7:18,7:18)=-eye(12);

L_x=[1 0 0; 0 0 0;0 0 0;0 0 0;0 0 1;0 1 0];
L_y=[0 0 0; 0 1 0;0 0 0;0 0 1;0 0 0;1 0 0];
L_z=[0 0 0; 0 0 0;0 0 1;0 1 0;1 0 0;0 0 0];

A_x(7:12,1:3)=C_hat*L_x;
A_y(7:12,1:3)=C_hat*L_y;
A_z(7:12,1:3)=C_hat*L_z;

A_x(13:18,4:6)=K_eq_til*C_tot*L_x;
A_y(13:18,4:6)=K_eq_til*C_tot*L_y;
A_z(13:18,4:6)=K_eq_til*C_tot*L_z;

R=M-1j*k_x*A_x-1j*k_y*A_y;
[V,D]=eig(A_z,R);
D=diag(D)*1j;
dof_S=[1 2 3 6 9 10 11 15];

[~,i_temp]=sort(real(D),'descend');

k_z=1./D(i_temp([1:4]));

SV=V(dof_S,i_temp([1:4 18:-1:15]));
SV(end,:)=-SV(end,:);


