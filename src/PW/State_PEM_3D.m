% State_PEM.m
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
% The State vector for a PEM material according to Dazel et al. JAP 2013
% S={\hat{\sigma}_{xz},u_z^s,u_z^t,\hat{\sigma}_{zz},p,u_x^s}
%
%
%

%%


function M=State_PEM_3D(k_x,k_y,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til)



k_z_1=sqrt(delta_1^2-k_x^2-k_y^2);
k_z_2=sqrt(delta_2^2-k_x^2-k_y^2);
k_z_3=sqrt(delta_3^2-k_x^2-k_y^2);


M(1:8,1)=[k_x;k_y; k_z_1; mu_1*k_z_1;-1j*(A_hat*delta_1^2+2*N*k_z_1^2);-2*1j*N*k_z_1*k_y;-2*1j*N*k_z_1*k_x;1j*delta_1^2*K_eq_til*mu_1];
M(1:8,5)=[k_x;k_y;-k_z_1;-mu_1*k_z_1;-1j*(A_hat*delta_1^2+2*N*k_z_1^2); 2*1j*N*k_z_1*k_y; 2*1j*N*k_z_1*k_x;1j*delta_1^2*K_eq_til*mu_1];
M(1:8,2)=[k_x;k_y; k_z_2; mu_2*k_z_2;-1j*(A_hat*delta_2^2+2*N*k_z_2^2);-2*1j*N*k_z_2*k_y;-2*1j*N*k_z_2*k_x;1j*delta_2^2*K_eq_til*mu_2];
M(1:8,6)=[k_x;k_y;-k_z_2;-mu_2*k_z_2;-1j*(A_hat*delta_2^2+2*N*k_z_2^2); 2*1j*N*k_z_2*k_y; 2*1j*N*k_z_2*k_x;1j*delta_2^2*K_eq_til*mu_2];
M(1:8,3)=[k_z_3;0;-k_x;-mu_3*k_x;2*1j*N*k_z_3*k_x;-1j*N*k_x*k_y;-1j*N*(k_z_3^2-k_x^2);0];
M(1:8,7)=[k_z_3;0; k_x; mu_3*k_x;2*1j*N*k_z_3*k_x; 1j*N*k_x*k_y; 1j*N*(k_z_3^2-k_x^2);0];
M(1:8,4)=[0;k_z_3;-k_y;-mu_3*k_y;2*1j*N*k_z_3*k_y;-1j*N*(k_z_3^2-k_y^2);-1j*N*k_x*k_y;0];
M(1:8,8)=[0;k_z_3; k_y; mu_3*k_y;2*1j*N*k_z_3*k_y; 1j*N*(k_z_3^2-k_y^2); 1j*N*k_x*k_y;0];






