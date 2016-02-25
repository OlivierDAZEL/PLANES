% interface_elas_PEM.m
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
eval(['Mat_elas_' num2str(medium_1-1000*floor(medium_1/1000))])
delta_P=omega*sqrt(rho_solide/(lambda_solide+2*mu_solide));
delta_s=omega*sqrt(rho_solide/mu_solide);
k_z_1=sqrt([delta_P delta_s].^2-k_x^2);
SV_1=State_elas_2D(k_x,delta_P,delta_s,lambda_solide,mu_solide);



eval(['Mat_porous_' num2str(medium_2-1000*floor(medium_2/1000))])
properties_eqf
properties_PEM
compute_Biot_waves
k_z_2=sqrt([delta_1 delta_2 delta_3].^2-k_x^2);
SV_2=State_PEM_2D(k_x,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til);




% sigma_xz=sigma_xz_hat
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_1(1))= SV_1(1,1)*exp(-1j*k_z_1(1)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(2))= SV_1(1,2)*exp(-1j*k_z_1(2)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(3))= SV_1(1,3);
Mat_PW(number_of_eq,dof_medium_1(4))= SV_1(1,4);
Mat_PW(number_of_eq,dof_medium_2(1))=-SV_2(1,1);
Mat_PW(number_of_eq,dof_medium_2(2))=-SV_2(1,2);
Mat_PW(number_of_eq,dof_medium_2(3))=-SV_2(1,3);
Mat_PW(number_of_eq,dof_medium_2(4))=-SV_2(1,4)*exp(-1j*k_z_2(1)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(5))=-SV_2(1,5)*exp(-1j*k_z_2(2)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(6))=-SV_2(1,6)*exp(-1j*k_z_2(3)*multilayer(i_interface+1).d);


% u_z=u_z_s
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_1(1))= SV_1(2,1)*exp(-1j*k_z_1(1)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(2))= SV_1(2,2)*exp(-1j*k_z_1(2)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(3))= SV_1(2,3);
Mat_PW(number_of_eq,dof_medium_1(4))= SV_1(2,4);
Mat_PW(number_of_eq,dof_medium_2(1))=-SV_2(2,1);
Mat_PW(number_of_eq,dof_medium_2(2))=-SV_2(2,2);
Mat_PW(number_of_eq,dof_medium_2(3))=-SV_2(2,3);
Mat_PW(number_of_eq,dof_medium_2(4))=-SV_2(2,4)*exp(-1j*k_z_2(1)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(5))=-SV_2(2,5)*exp(-1j*k_z_2(2)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(6))=-SV_2(2,6)*exp(-1j*k_z_2(3)*multilayer(i_interface+1).d);

% u_z=u_z_t
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_1(1))= SV_1(2,1)*exp(-1j*k_z_1(1)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(2))= SV_1(2,2)*exp(-1j*k_z_1(2)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(3))= SV_1(2,3);
Mat_PW(number_of_eq,dof_medium_1(4))= SV_1(2,4);
Mat_PW(number_of_eq,dof_medium_2(1))=-SV_2(3,1);
Mat_PW(number_of_eq,dof_medium_2(2))=-SV_2(3,2);
Mat_PW(number_of_eq,dof_medium_2(3))=-SV_2(3,3);
Mat_PW(number_of_eq,dof_medium_2(4))=-SV_2(3,4)*exp(-1j*k_z_2(1)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(5))=-SV_2(3,5)*exp(-1j*k_z_2(2)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(6))=-SV_2(3,6)*exp(-1j*k_z_2(3)*multilayer(i_interface+1).d);


% sigma_zz=sigma_zz_t-p
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_1(1))= SV_1(3,1)*exp(-1j*k_z_1(1)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(2))= SV_1(3,2)*exp(-1j*k_z_1(2)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(3))= SV_1(3,3);
Mat_PW(number_of_eq,dof_medium_1(4))= SV_1(3,4);
Mat_PW(number_of_eq,dof_medium_2(1))=-(SV_2(4,1)-SV_2(5,1));
Mat_PW(number_of_eq,dof_medium_2(2))=-(SV_2(4,2)-SV_2(5,2));
Mat_PW(number_of_eq,dof_medium_2(3))=-(SV_2(4,3)-SV_2(5,3));
Mat_PW(number_of_eq,dof_medium_2(4))=-(SV_2(4,4)-SV_2(5,4))*exp(-1j*k_z_2(1)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(5))=-(SV_2(4,5)-SV_2(5,5))*exp(-1j*k_z_2(2)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(6))=-(SV_2(4,6)-SV_2(5,6))*exp(-1j*k_z_2(3)*multilayer(i_interface+1).d);


% u_x=u_x_s
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_1(1))= SV_1(4,1)*exp(-1j*k_z_1(1)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(2))= SV_1(4,2)*exp(-1j*k_z_1(2)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(3))= SV_1(4,3);
Mat_PW(number_of_eq,dof_medium_1(4))= SV_1(4,4);
Mat_PW(number_of_eq,dof_medium_2(1))=-SV_2(6,1);
Mat_PW(number_of_eq,dof_medium_2(2))=-SV_2(6,2);
Mat_PW(number_of_eq,dof_medium_2(3))=-SV_2(6,3);
Mat_PW(number_of_eq,dof_medium_2(4))=-SV_2(6,4)*exp(-1j*k_z_2(1)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(5))=-SV_2(6,5)*exp(-1j*k_z_2(2)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(6))=-SV_2(6,6)*exp(-1j*k_z_2(3)*multilayer(i_interface+1).d);

