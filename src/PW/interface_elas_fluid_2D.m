% interface_elas_fluid.m
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

if medium_2==0
    k_z_2=sqrt(k_air^2-k_x^2);
    SV_2=State_fluid_2D(k_x,k_z_2,air.K);
else
    eval(['Mat_porous_' num2str(medium_2-1000*floor(medium_2/1000))])
    typ_mat=floor(medium_2/1000);
    properties_eqf
    k_z_2=sqrt((omega*sqrt(rho_eq_til/K_eq_til))^2-k_x^2);
    SV_2=State_fluid_2D(k_x,k_z_2,K_eq_til);
end


% sigma_xz=0
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_1(1))= SV_1(1,1)*exp(-1j*k_z_1(1)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(2))= SV_1(1,2)*exp(-1j*k_z_1(2)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(3))= SV_1(1,3);
Mat_PW(number_of_eq,dof_medium_1(4))= SV_1(1,4);

% u_z=u_z
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_1(1))= SV_1(2,1)*exp(-1j*k_z_1(1)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(2))= SV_1(2,2)*exp(-1j*k_z_1(2)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(3))= SV_1(2,3);
Mat_PW(number_of_eq,dof_medium_1(4))= SV_1(2,4);
Mat_PW(number_of_eq,dof_medium_2(1))=-SV_2(1,1);
Mat_PW(number_of_eq,dof_medium_2(2))=-SV_2(1,2)*exp(-1j*k_z_2*multilayer(i_interface+1).d);


% sigma_zz=-p
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_1(1))= SV_1(3,1)*exp(-1j*k_z_1(1)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(2))= SV_1(3,2)*exp(-1j*k_z_1(2)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(3))= SV_1(3,3);
Mat_PW(number_of_eq,dof_medium_1(4))= SV_1(3,4);
Mat_PW(number_of_eq,dof_medium_2(1))= SV_2(2,1);
Mat_PW(number_of_eq,dof_medium_2(2))= SV_2(2,2)*exp(-1j*k_z_2*multilayer(i_interface+1).d);

