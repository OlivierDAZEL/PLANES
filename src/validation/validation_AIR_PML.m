% validation_AIR_PML.m
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




tau_x=exp(j*pi/4);
%tau_x=1

M_an=[];
F_an=[];


k_air=omega/air.c;

x=-d_PML;
M_an(1,1)=sin(k_air*tau_x*x);
M_an(1,2)=-sin(k_air*x);
M_an(1,3)=-cos(k_air*x);

M_an(2,1)=air.K*k_air*cos(k_air*tau_x*x);
M_an(2,2)=-air.K*k_air*cos(k_air*x);
M_an(2,3)= air.K*k_air*sin(k_air*x);

x=-d_PML-d_air;
M_an(3,2)=sin(k_air*x);
M_an(3,3)=cos(k_air*x);

%F_an(2,1)=-1;
F_an(3,1)=-1;

x_poreux=linspace(-(d_PML),0,100);
x_air=linspace(-(d_air+d_PML),-d_PML,100);

X_an=M_an\F_an;

A=X_an(1);
B=X_an(2);
C=X_an(3);

u_x_poreux=A*sin(k_air*tau_x*x_poreux);
v_z_poreux=0*u_x_poreux;
v_x_poreux=j*omega*u_x_poreux;
sigma_p=air.K*k_air*A*cos(k_air*tau_x*x_poreux);

u_x_air=B*sin(k_air*x_air)+C*cos(k_air*x_air);
v_z_air=0*u_x_air;
v_x_air=j*omega*u_x_air;
sigma_a=air.K*k_air*(B*cos(k_air*x_air)-C*sin(k_air*x_air));




% figure(1000)
% subplot(221)
hold on
plot(x_air+d_air+d_PML,abs(sigma_a))
plot(x_poreux+d_air+d_PML,abs(sigma_p),'r')
% subplot(223)
% hold on
% plot(x_air,abs(sigma_a))
% plot(x_poreux,abs(sigma_p),'r')
% 
% figure(1001)
% subplot(221)
% hold on
% plot(x_air,angle(v_x_air))
% plot(x_poreux,angle(v_x_poreux),'r')
% subplot(223)
% hold on
% plot(x_air,angle(sigma_a))
% plot(x_poreux,angle(sigma_p),'r')

