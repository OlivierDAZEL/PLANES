% plot_solution_TR6.m
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




for ie=1:nb_elements
     coord_e=nodes(elements(ie,:),:);
     sol_e=sol(elements(ie,:));
     xdata=[coord_e( [1 2 6],1) coord_e( [2 3 4],1) coord_e( [2 4 6],1) coord_e( [4 5 6],1)];
     ydata=[coord_e( [1 2 6],2) coord_e( [2 3 4],2) coord_e( [2 4 6],2) coord_e( [4 5 6],2)];
     zdata=abs([sol_e( [1 2 6]) sol_e( [2 3 4]) sol_e( [2 4 6]) sol_e( [4 5 6])]);
     patch(xdata,ydata,zdata,'EdgeColor','none')
     line([coord_e([1 3 5 1],1)],[coord_e([1 3 5 1],2)],'Color','k')
end
     

