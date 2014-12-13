% TR6_elas.m
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


function [vm,vk0,vk1] = TR6_elas(coord_e)% Gauss Points Data%------------------------------npg=6;a=0.445948490915965;b=0.091576213509771;P1=0.111690794839005;P2=0.054975871827661;ksig=[  a, a;    1-2*a, a;    a,   1-2*a;    b, b;    1-2*b, b;    b, 1-2*b];w_g=[P1,P1,P1,P2,P2,P2];% Matrices initialization%-----------------------vk0 = zeros(12,12);vk1 = zeros(12,12);vm=zeros(12,12);D0=[1 1 0;1 1 0; 0 0 0];D1=[0 -2 0;-2 0 0; 0 0 1];% Loop on Gauss Points%------------------for ipg=1:npg,    ksi=ksig(ipg,1);    eta=ksig(ipg,2);    lambda=1-ksi-eta;        Phi =  [ -lambda*(1-2*lambda)  4*ksi*lambda  -ksi*(1-2*ksi) 4*ksi*eta -eta*(1-2*eta) 4*eta*lambda];    dPhi  = [ 1-4*lambda 4*(lambda-ksi) -1+4*ksi 4*eta 0 -4*eta; ...             1-4*lambda -4*ksi 0 4*ksi -1+4*eta 4*(lambda-eta)];    J = dPhi*coord_e';        weight = w_g(ipg) * det(J);        vB =J\dPhi;    eps=[vB(1,1) 0 vB(1,2) 0 vB(1,3) 0 vB(1,4) 0 vB(1,5) 0 vB(1,6) 0;           0 vB(2,1) 0 vB(2,2) 0 vB(2,3) 0 vB(2,4) 0 vB(2,5) 0 vB(2,6);           vB(2,1) vB(1,1) vB(2,2) vB(1,2) vB(2,3) vB(1,3) vB(2,4) vB(1,4) vB(2,5) vB(1,5) vB(2,6) vB(1,6)];    Phi_u=[Phi(1) 0 Phi(2) 0 Phi(3) 0 Phi(4) 0 Phi(5) 0 Phi(6) 0;        0 Phi(1) 0 Phi(2) 0 Phi(3) 0 Phi(4) 0 Phi(5) 0 Phi(6)];        vk0= vk0 +  eps'*D0*eps *weight;    vk1= vk1 +  eps'*D1*eps *weight;    vm = vm +  Phi_u'*Phi_u * weight;end;   % fin de boucle sur les points de gauss.