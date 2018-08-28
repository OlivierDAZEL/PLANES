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
if mat_typ==2 || mat_typ==3
    eval(['Mat_porous_' num2str(media(1).mat-1000*mat_typ)])
end

switch (mat_typ)
case {0}
    K_fluid = air.K;
    rho_fluid = air.rho;
case {2}
    properties_JCA;
    K_fluid = K_eq_til;
    rho_fluid = rho_eq_til;
case {3}
    properties_JCA
    eqf2limp
    K_fluid = K_eq_til;
    rho_fluid = rho_eq_til;
end

d=signe*media(1).d;

A=zeros(2,2);
A(1,2) = -1/K_fluid + k_x^2/(rho_fluid*omega^2);
A(2,1) = rho_fluid*omega^2;

MM=expm(A*d);

Omega_plus=MM*Omega_moins;
