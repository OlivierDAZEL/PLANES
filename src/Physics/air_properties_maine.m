% air_properties_maine.m
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






%% Constants
air.T=20;                                     %[?C]
air.P=0.10132e6;                              %[Pa]
air.gamma=1.400;                                %[1]
air.mmol=.29e-1;                                %[kg.mol^-1]
air.mu=.184e-4; % Dynamic viscosity           %(kg/m/s)
air.NPR=0.710;
air.lambda=0.0262;                              %[W.m^-1.K^-1]

%% Air properties
air.rho=1.204;                                %[kg.m^-3]
air.c=343.25946;                              %[m.s^-1]

air.K=air.c^2*air.rho;
air.Z=air.rho*air.c;

air.C_p=1006;
air.C_v=air.C_p/air.gamma;

air.nu=air.mu/air.rho;
air.nu_prime=air.nu/air.NPR;