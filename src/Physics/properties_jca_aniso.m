% properties_jca.m
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

Q=rotate_u3(angle_x,angle_y,angle_z);

sig_tensor=Q*sig_tensor_0*Q';
LCV_tensor=Q*LCV_tensor_0*Q';
alpha_tensor=Q*alpha_tensor_0*Q';



%%  Johnson et al model for rho_eq_til


% omega_0=sig_tensor*phi*inv(air.rho*alpha_tensor);
% omega_infty=(sig_tensor.*phi.*LCV_tensor).^2*inv(4*air.mu*air.rho.*alpha_tensor.^2);
% F_JKD=sqrt(eye(3)+1j*omega*inv(omega_infty));
% rho_eq_til=(air.rho*alpha_tensor/phi)*(eye(3)+(omega_0/(1j*omega))*F_JKD);
% alpha_til_tensor=phi*rho_eq_til./air.rho;


q_0=air.mu*inv(sig_tensor);
alpha_til=alpha_tensor+((air.mu*phi)/(1i*omega*air.rho*air.mu))*sig_tensor*sqrtm(eye(3)+((1j*4*omega*air.rho)/(air.mu*phi^2))*alpha_tensor^2*q_0^2*inv(LCV_tensor).^2);
rho_eq_til=(air.rho/phi)*alpha_til;




%%  Champoux-Allard model for K_eq_til

omega_prime_infty=(16*air.nu_prime)/(LCT^2);
F_prime_CA=sqrt(1+1j*omega./omega_prime_infty);
alpha_prime_til=1+omega_prime_infty.*F_prime_CA./(2*1j*omega);
K_eq_til=(air.gamma*air.P./phi)./(air.gamma-(air.gamma-1)./alpha_prime_til);