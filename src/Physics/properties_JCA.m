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




%%  Johnson et al model for rho_eq_til

omega_0=sig*phi/(air.rho*alpha);
omega_infty=(sig*phi*LCV)^2/(4*air.mu*air.rho*alpha^2);
F_JKD=sqrt(1+1j*omega/omega_infty);
rho_eq_til=(air.rho*alpha/phi)*(1+(omega_0/(1j*omega))*F_JKD);
alpha_til=phi*rho_eq_til./air.rho;


%%  Champoux-Allard model for K_eq_til

omega_prime_infty=(16*air.nu_prime)/(LCT^2);
F_prime_CA=sqrt(1+1j*omega./omega_prime_infty);
alpha_prime_til=1+omega_prime_infty.*F_prime_CA./(2*1j*omega);
K_eq_til=(air.gamma*air.P./phi)./(air.gamma-(air.gamma-1)./alpha_prime_til);



