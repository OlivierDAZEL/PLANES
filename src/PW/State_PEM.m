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
%%


function M=State_PEM(k_x,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til)



beta_1=sqrt(delta_1^2-k_x^2);
beta_2=sqrt(delta_2^2-k_x^2);
beta_3=sqrt(delta_3^2-k_x^2);

alpha_1=-1i*A_hat*delta_1^2-1i*2*N*beta_1^2;
alpha_2=-1i*A_hat*delta_2^2-1i*2*N*beta_2^2;
alpha_3= 2*1i*N*beta_3*k_x;


M(1:6,1)=[-2*1i*N*beta_1*k_x; beta_1; mu_1*beta_1;alpha_1;1i*delta_1^2*K_eq_til*mu_1;k_x];
M(1:6,4)=[ 2*1i*N*beta_1*k_x;-beta_1;-mu_1*beta_1;alpha_1;1i*delta_1^2*K_eq_til*mu_1;k_x];

M(1:6,2)=[-2*1i*N*beta_2*k_x; beta_2; mu_2*beta_2;alpha_2;1i*delta_2^2*K_eq_til*mu_2;k_x];
M(1:6,5)=[ 2*1i*N*beta_2*k_x;-beta_2;-mu_2*beta_2;alpha_2;1i*delta_2^2*K_eq_til*mu_2;k_x];

M(1:6,3)=[1i*N*(beta_3^2-k_x^2);k_x;mu_3*k_x;-alpha_3;0;-beta_3];
M(1:6,6)=[1i*N*(beta_3^2-k_x^2);k_x;mu_3*k_x; alpha_3;0; beta_3];