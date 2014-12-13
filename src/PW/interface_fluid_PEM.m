% interface_fluid_PEM.m
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
if medium_1==0
    k_z_1=sqrt(k_air^2-k_x^2);
    SV_1=State_fluid(k_x,k_z_1,air.K);
else
    eval(['Mat_fluid_' num2str(medium_1-1000*floor(medium_1/1000))])
    properties_jca
    k_z_1=sqrt((omega*sqrt(rho_eq_til/K_eq_til))^2-k_x^2);
    SV_1=State_fluid(k_x,k_z_1,K_eq_til);
end
eval(['Mat_PEM_' num2str(medium_2-1000*floor(medium_2/1000))])
properties_jca
properties_PEM
compute_Biot_waves
k_z_2=sqrt([delta_1 delta_2 delta_3].^2-k_x^2);
SV_2=State_PEM(k_x,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til);



% Continuity of normal displacement
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_1(1))= SV_1(1,1)*exp(-1j*k_z_1*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(2))= SV_1(1,2);
Mat_PW(number_of_eq,dof_medium_2(1))=-SV_2(3,1);
Mat_PW(number_of_eq,dof_medium_2(2))=-SV_2(3,2);
Mat_PW(number_of_eq,dof_medium_2(3))=-SV_2(3,3);
Mat_PW(number_of_eq,dof_medium_2(4))=-SV_2(3,4)*exp(-1j*k_z_2(1)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(5))=-SV_2(3,5)*exp(-1j*k_z_2(2)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(6))=-SV_2(3,6)*exp(-1j*k_z_2(3)*multilayer(i_interface+1).d);

% Continuity of pressure
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_1(1))= SV_1(2,1)*exp(-1j*k_z_1*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(2))= SV_1(2,2);
Mat_PW(number_of_eq,dof_medium_2(1))=-SV_2(5,1);
Mat_PW(number_of_eq,dof_medium_2(2))=-SV_2(5,2);
Mat_PW(number_of_eq,dof_medium_2(3))=-SV_2(5,3);
Mat_PW(number_of_eq,dof_medium_2(4))=-SV_2(5,4)*exp(-1j*k_z_2(1)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(5))=-SV_2(5,5)*exp(-1j*k_z_2(2)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(6))=-SV_2(5,6)*exp(-1j*k_z_2(3)*multilayer(i_interface+1).d);


% Nullity of sigma_xx
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_2(1))=SV_2(1,1);
Mat_PW(number_of_eq,dof_medium_2(2))=SV_2(1,2);
Mat_PW(number_of_eq,dof_medium_2(3))=SV_2(1,3);
Mat_PW(number_of_eq,dof_medium_2(4))=SV_2(1,4)*exp(-1j*k_z_2(1)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(5))=SV_2(1,5)*exp(-1j*k_z_2(2)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(6))=SV_2(1,6)*exp(-1j*k_z_2(3)*multilayer(i_interface+1).d);


% Nullity of sigma_zz
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_2(1))=SV_2(4,1);
Mat_PW(number_of_eq,dof_medium_2(2))=SV_2(4,2);
Mat_PW(number_of_eq,dof_medium_2(3))=SV_2(4,3);
Mat_PW(number_of_eq,dof_medium_2(4))=SV_2(4,4)*exp(-1j*k_z_2(1)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(5))=SV_2(4,5)*exp(-1j*k_z_2(2)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(6))=SV_2(4,6)*exp(-1j*k_z_2(3)*multilayer(i_interface+1).d);


