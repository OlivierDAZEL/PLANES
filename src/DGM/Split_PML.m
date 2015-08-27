% Split_PML.m
%
% Copyright (C) 2015 < Olivier DAZEL <olivier.dazel@univ-lemans.fr> >
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
%                 https://github.com/OlivierDAZEL/PLANES
% or find more details on Olivier's webpage
%                  http://perso.univ-lemans.fr/~odazel/
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
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [F_plus,F_moins]=Split_PML(nx,ny,Z_e,tau_x,tau_y)

ne=sqrt((nx/tau_x)^2+(ny/tau_y)^2);


F_plus(1,1)=-Z_e*(nx^2/2)/(ne*tau_x^2);
F_plus(1,2)=-Z_e*(nx*ny/2)/(ne*tau_x*tau_y);
F_plus(1,3)=(nx/2)/tau_x;
F_plus(2,1)=-Z_e*(nx*ny/2)/(ne*tau_x*tau_y);
F_plus(2,2)=-Z_e*(ny^2/2)/(ne*tau_y^2);
F_plus(2,3)=(ny/2)/tau_y;
F_plus(3,1)=(nx/2)/tau_x;
F_plus(3,2)=(ny/2)/tau_y;
F_plus(3,3)=-ne*(0.5)/Z_e;


F_moins(1,1)=Z_e*(nx^2/2)/(ne*tau_x^2);
F_moins(1,2)=Z_e*(nx*ny/2)/(ne*tau_x*tau_y);
F_moins(1,3)=(nx/2)/tau_x;
F_moins(2,1)=Z_e*(nx*ny/2)/(ne*tau_x*tau_y);
F_moins(2,2)=Z_e*(ny^2/2)/(ne*tau_y^2);
F_moins(2,3)=(ny/2)/tau_y;
F_moins(3,1)=(nx/2)/tau_x;
F_moins(3,2)=(ny/2)/tau_y;
F_moins(3,3)=ne*(0.5)/Z_e;






