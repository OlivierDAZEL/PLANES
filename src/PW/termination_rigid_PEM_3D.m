% termination_rigid_PEM.m
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


% u_x^s=0
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_2(1))=SV_2(1,1)*exp(-1j*k_z_2(1)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(2))=SV_2(1,2)*exp(-1j*k_z_2(2)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(3))=SV_2(1,3)*exp(-1j*k_z_2(3)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(4))=SV_2(1,4)*exp(-1j*k_z_2(4)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(5))=SV_2(1,5);
Mat_PW(number_of_eq,dof_medium_2(6))=SV_2(1,6);
Mat_PW(number_of_eq,dof_medium_2(7))=SV_2(1,7);
Mat_PW(number_of_eq,dof_medium_2(8))=SV_2(1,8);

% u_x^s=0
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_2(1))=SV_2(2,1)*exp(-1j*k_z_2(1)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(2))=SV_2(2,2)*exp(-1j*k_z_2(2)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(3))=SV_2(2,3)*exp(-1j*k_z_2(3)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(4))=SV_2(2,4)*exp(-1j*k_z_2(4)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(5))=SV_2(2,5);
Mat_PW(number_of_eq,dof_medium_2(6))=SV_2(2,6);
Mat_PW(number_of_eq,dof_medium_2(7))=SV_2(2,7);
Mat_PW(number_of_eq,dof_medium_2(8))=SV_2(2,8);

% u_x^s=0
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_2(1))=SV_2(3,1)*exp(-1j*k_z_2(1)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(2))=SV_2(3,2)*exp(-1j*k_z_2(2)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(3))=SV_2(3,3)*exp(-1j*k_z_2(3)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(4))=SV_2(3,4)*exp(-1j*k_z_2(4)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(5))=SV_2(3,5);
Mat_PW(number_of_eq,dof_medium_2(6))=SV_2(3,6);
Mat_PW(number_of_eq,dof_medium_2(7))=SV_2(3,7);
Mat_PW(number_of_eq,dof_medium_2(8))=SV_2(3,8);

% u_x^s=0
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_2(1))=SV_2(4,1)*exp(-1j*k_z_2(1)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(2))=SV_2(4,2)*exp(-1j*k_z_2(2)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(3))=SV_2(4,3)*exp(-1j*k_z_2(3)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(4))=SV_2(4,4)*exp(-1j*k_z_2(4)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(5))=SV_2(4,5);
Mat_PW(number_of_eq,dof_medium_2(6))=SV_2(4,6);
Mat_PW(number_of_eq,dof_medium_2(7))=SV_2(4,7);
Mat_PW(number_of_eq,dof_medium_2(8))=SV_2(4,8);
