% Mat_PEM_3.m
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

porous_model.eqf='JCA_aniso';
porous_model.frame='structural';
porous_model.aniso='yes';

angle_rot=[pi/3;4*pi/9;pi/4];
Q=rotate_u3(angle_rot);

phi=0.95;
LCT=4.500E-05;
rho_1=126.000;

young=694400E+00;
nu=0.24000E+00;
eta=0.05;


sig=Q*diag([10000;20000;40000])*Q';
alpha=Q*diag([1.1;1.1;1.1])*Q';
LCV=Q*diag([1.50E-05;1.50E-05;1.50E-05])*Q';



N=(young)/(2*(1+nu));
A_hat=(young*nu)/((1+nu)*(1-2*nu));
F_hat=1*A_hat;
C_hat_0= [A_hat+2*N A_hat F_hat 0 0 0;A_hat A_hat+2*N F_hat 0 0 0;F_hat F_hat A_hat+2*N 0 0 0; 0 0 0 N 0 0;0 0 0 0 N 0; 0 0 0 0 0 N ];



alpha_hat=0.33348;
beta_hat =812.69e3;
b_hat=0.29620;
