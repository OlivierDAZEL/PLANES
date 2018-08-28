% transfer_unknowns.m
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

function [Omega_plus]=transfer_unknowns(k_x,omega,Omega_moins,signe,media,air)

mat_typ = floor(media(1).mat/1000);
eval(['Mat_porous_' num2str(media(1).mat-1000*mat_typ)])
properties_JCA
properties_PEM
d=signe*media(1).d;

A=zeros(6,6);
A(1,4) = 1j*k_x*A_hat/P_hat;
A(1,5) = 1j*gamma_til*k_x;
A(1,6) = -(A_hat^2-P_hat^2)/P_hat*k_x^2-rho_til*omega^2;

A(2,4) = 1/P_hat;
A(2,6) = 1j*k_x*A_hat/P_hat;

A(3,5) = -1/K_eq_til + k_x^2/(rho_eq_til*omega^2);
A(3,6) = -1j*k_x*gamma_til;

A(4,1) = 1j*k_x;
A(4,2) = -rho_s_til*omega^2;
A(4,3) = -rho_eq_til*gamma_til*omega^2;

A(5,2) = rho_eq_til*gamma_til*omega^2;
A(5,3) = rho_eq_til*omega^2;

A(6,1) = 1/N;
A(6,2) = 1j*k_x;

MM=expm(A*d);

Omega_plus=MM*Omega_moins;
