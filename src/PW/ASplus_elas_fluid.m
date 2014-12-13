% ASplus_elas_fluid.m
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


% Initialization of the State vectors in media 1 and 2

eval(['Mat_elas_' num2str(multilayer(end).mat-1000*floor(multilayer(end).mat/1000))])
delta_P=omega*sqrt(rho_solide/(lambda_solide+2*mu_solide));
delta_s=omega*sqrt(rho_solide/mu_solide);
k_z_2=sqrt([delta_P delta_s].^2-k_x^2);
dof_medium_2=nS_minus+nb_amplitudes+1-[4:-1:1];



% Continuity of normal displacement
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,nS_minus+nb_amplitudes+1)= -1;
Mat_PW(number_of_eq,dof_medium_2(1))=SV_2(2,1)*exp(-1j*k_z_2(1)*multilayer(end).d);
Mat_PW(number_of_eq,dof_medium_2(2))=SV_2(2,2)*exp(-1j*k_z_2(2)*multilayer(end).d);
Mat_PW(number_of_eq,dof_medium_2(3))=SV_2(2,3);
Mat_PW(number_of_eq,dof_medium_2(4))=SV_2(2,4);
% Continuity of the pressure
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,nS_minus+nb_amplitudes+2)= 1;
Mat_PW(number_of_eq,dof_medium_2(1))=SV_2(3,1)*exp(-1j*k_z_2(1)*multilayer(end).d);
Mat_PW(number_of_eq,dof_medium_2(2))=SV_2(3,2)*exp(-1j*k_z_2(2)*multilayer(end).d);
Mat_PW(number_of_eq,dof_medium_2(3))=SV_2(3,3);
Mat_PW(number_of_eq,dof_medium_2(4))=SV_2(3,4);

% Nullity of sigma_xz
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_2(1))=SV_2(1,1)*exp(-1j*k_z_2(1)*multilayer(end).d);
Mat_PW(number_of_eq,dof_medium_2(2))=SV_2(1,2)*exp(-1j*k_z_2(2)*multilayer(end).d);
Mat_PW(number_of_eq,dof_medium_2(3))=SV_2(1,3);
Mat_PW(number_of_eq,dof_medium_2(4))=SV_2(1,4);
