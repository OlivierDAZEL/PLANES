
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
% S={{\sigma}_{xz},u_z,{\sigma}_{zz},p,u_x}
%
%
%
%%


function M=State_elas(k_x,delta_P,delta_s,lambda,mu)

% Generate the State vector
% SV= [\sigma_xz u_z \sigma_zz u_x]

beta_P=sqrt(delta_P^2-k_x^2);
beta_s=sqrt(delta_s^2-k_x^2);

alpha_P=-1j*lambda*delta_P^2-1j*2*mu*beta_P^2;
alpha_s= 2*1j*mu*beta_s*k_x;

M(1:4,1)=[-2*1j*mu*beta_P*k_x; beta_P;alpha_P;k_x];
M(1:4,3)=[ 2*1j*mu*beta_P*k_x;-beta_P;alpha_P;k_x];

M(1:4,2)=[1j*mu*(beta_s^2-k_x^2);k_x;-alpha_s;-beta_s];
M(1:4,4)=[1j*mu*(beta_s^2-k_x^2);k_x; alpha_s; beta_s];
