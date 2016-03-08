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
porous_model.frame='anelastic';
porous_model.aniso='yes';



angle_rot=[0;3*pi/4;0];
Q=rotate_u3(angle_rot);

rho_1=9.200;
phi=0.99;
LCV=240E-06;
LCT=490E-06;
alpha=1.2;


eigvec=[   0.97000   0.19000   0.12000;
  -0.14000  -0.92000   0.36000;
   0.18000   0.34000   0.92000];
eigvalues=[9800 0 0; 0 10500 0; 0 0 11400];
 
sig_i=diag(diag(inv(eigvec)*(eigvalues*eigvec)));



sig=Q*sig_i*Q';
alpha=Q*diag([1.2;1.2;1.2])*Q';
LCV=Q*diag([240E-06;240E-06;240E-06])*Q';


C_hat_0=1e5*[
   7.719353596027278   3.425199931688752  -0.022642853370394                   0                   0                   0
   3.425199931688752   4.278234531826723   1.184528821838865                   0                   0                   0
  -0.022642853370394   1.184528821838865   2.215513353247722                   0                   0                   0
                   0                   0                   0   1.036365276907575                   0                   0
                   0                   0                   0                   0   1.236840201799869                   0
                   0                   0                   0                   0                   0   1.012310180008057];


alpha_hat=0.33348;
beta_hat =832.16e3;
b_hat=0.29620;
