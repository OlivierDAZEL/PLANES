% State_fluid_3D.m
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


function [k_z,SV]=State_general_3D(k_x,k_y,omega,air)


Mat_porous_3
properties_eqf
properties_PEM
compute_Biot_waves
% k_z_1=sqrt([delta_1 delta_2 delta_3 delta_3].^2-k_x^2);
% SV_1=State_PEM_3D(k_x,k_y,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til);

M=zeros(13,13);
A_x=zeros(13,13);
A_y=zeros(13,13);
A_z=zeros(13,13);


M(1:3,1:3)=-omega^2*rho_s_til*eye(3);
M(1:3,4:6)=-omega^2*gamma_til*rho_eq_til*eye(3);

A_x(1,7)=-1;
A_y(1,12)=-1;
A_z(1,11)=-1;

A_x(2,10)=-1;
A_y(2,8)=-1;
A_z(2,12)=-1;

A_x(3,11)=-1;
A_y(3,10)=-1;
A_z(3,9)=-1;



M(4:6,1:3)=-omega^2*gamma_til*rho_eq_til*eye(3);
M(4:6,4:6)=-omega^2*rho_eq_til*eye(3);

A_x(4,13)=1;
A_y(5,13)=1;
A_z(6,13)=1;

M(7:12,7:12)=-eye(6);


L_x=[1 0 0; 0 0 0;0 0 0;0 0 0;0 0 1;0 0 1];
L_y=[0 0 0; 0 1 0;0 0 0;0 0 1;0 0 0;1 0 0];
L_z=[0 0 0; 0 0 0;0 0 1;0 1 0;1 0 0;0 0 0];

A_x(7:12,1:3)=C_hat*L_x;
A_y(7:12,1:3)=C_hat*L_y;
A_z(7:12,1:3)=C_hat*L_z;



M(13,13)=1;
A_x(13,4)=K_eq_til;
A_y(13,5)=K_eq_til;
A_z(13,6)=K_eq_til;


R=M-1j*k_x*A_x-1j*k_y*A_y;
[V,D]=eig(A_z,R);
D=diag(D)*1j;
dof_S=[1 2 3 6 9 10 11 13];
dof_S_prime=[1 2];
length_S=8;
length_S_prime=5;
[~,i_temp]=sort(real(D),'descend');
k_z=1./D(i_temp([1:4]));
SV=V(dof_S,i_temp([1:4 13:-1:10]));

