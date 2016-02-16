% State_elas.m
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
% The State vector for a PEM material according to Dazel et al. JAP 2013
% S={u_x,u_y,u_z,\sigma_{zz},\sigma_{yz},\sigma_{xz}}
%
%
%
%%


function M=State_elas_3D(k_x,k_y,delta_P,delta_s,lambda,mu)

k_z_P=sqrt(delta_P^2-k_x^2-k_y^2);
k_z_S=sqrt(delta_s^2-k_x^2-k_y^2);


%  |k_x     exp(-j(k_x x + k_y y +k_z z))
%  |k_y
%  |k_z_P
%

M(1:6,1)=[k_x;k_y; k_z_P;-1j*(lambda*delta_P^2+2*mu*k_z_P^2);-2*1j*mu*k_z_P*k_y;-2*1j*mu*k_z_P*k_x];
M(1:6,4)=[k_x;k_y;-k_z_P;-1j*(lambda*delta_P^2+2*mu*k_z_P^2); 2*1j*mu*k_z_P*k_y; 2*1j*mu*k_z_P*k_x];

M(1:6,2)=[k_z_S;0;-k_x;2*1j*mu*k_z_S*k_x;-1j*mu*k_x*k_y;-1j*mu*(k_z_S^2-k_x^2)];
M(1:6,5)=[k_z_S;0; k_x;2*1j*mu*k_z_S*k_x; 1j*mu*k_x*k_y; 1j*mu*(k_z_S^2-k_x^2)];
 
M(1:6,3)=[0;k_z_S;-k_y;2*1j*mu*k_z_S*k_y;-1j*mu*(k_z_S^2-k_y^2);-1j*mu*k_x*k_y];
M(1:6,6)=[0;k_z_S; k_y;2*1j*mu*k_z_S*k_y; 1j*mu*(k_z_S^2-k_y^2); 1j*mu*k_x*k_y];


