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

function [Omega_plus]=transfer_unknowns(k_x,omega,Omega_moins,signe,media)

eval(['Mat_elas_' num2str(media(1).mat-1000*floor(media(1).mat/1000))])


lambda=lambda_solide;
mu=mu_solide;
rho=rho_solide;
d=signe*media(1).d;

P=lambda+2*mu;

A=zeros(4,4);
A(1,3)=1j*k_x*lambda/P;
A(1,4)=(-(lambda^2-P^2)/P)*k_x^2-rho*omega^2;

A(2,3)=1/P;
A(2,4)=1j*k_x*lambda/P;

A(3,1)=1j*k_x;
A(3,2)=-rho*omega^2;

A(4,1)=1/mu;
A(4,2)=1j*k_x;

MM=expm(A*d);

Omega_plus=MM*Omega_moins;
