% integrate_polynom_2D_edge.m
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
%%



function f = integrate_polynom_2D_edge(P1,base_e1,P2,base_e2,a,b,Gauss_points)

delta=(b-a);
points=((a+b)/2)*ones(1,Gauss_points.nb)+delta*Gauss_points.xi/2;

points_e1=points-base_e1'*ones(1,Gauss_points.nb);
points_e2=points-base_e2'*ones(1,Gauss_points.nb);
temp1=evaluate_polynom_2D_vectorial(P1,points_e1(1,:),points_e1(2,:));
temp2=evaluate_polynom_2D_vectorial(P2,points_e2(1,:),points_e2(2,:));
temp=sum(temp1.*temp2.*Gauss_points.w);
f=temp*norm(delta)/2;

end



