% compute_Biot_waves.m
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


delta_eq=omega.*sqrt(rho_eq_til./K_eq_til);
delta_s_1=omega.*sqrt(rho_til./P_hat);
delta_s_2=omega.*sqrt(rho_s_til./P_hat);


Psi=((delta_s_2.^2+delta_eq.^2).^2-4*delta_eq.^2.*delta_s_1.^2);
sdelta_total=sqrt(Psi);
delta_1=sqrt(0.5*(delta_s_2.^2+delta_eq.^2+sdelta_total));
delta_2=sqrt(0.5*(delta_s_2.^2+delta_eq.^2-sdelta_total));



mu_1=gamma_til.*delta_eq.^2./(delta_1.^2-delta_eq.^2);
mu_2=gamma_til.*delta_eq.^2./(delta_2.^2-delta_eq.^2);

delta_3=omega.*sqrt(rho_til/N);
mu_3=-gamma_til;

