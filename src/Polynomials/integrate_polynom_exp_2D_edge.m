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



function f = integrate_polynom_exp_2D_edge(P,base_polynom,jk,center_wave,a,b,Gauss_points)


delta_gauss=(b-a);

points=((a+b)/2)*ones(1,Gauss_points.nb)+delta_gauss*Gauss_points.xi/2;

points_poly=points-base_polynom'*ones(1,Gauss_points.nb);


temp=evaluate_polynom_2D_vectorial(P,points_poly(1,:),points_poly(2,:));

temp=sum(temp.*exp(jk.'*(points-center_wave*ones(1,Gauss_points.nb))).*Gauss_points.w);
f=temp*norm(delta_gauss)/2;


end



