% integrate_exp_element.m
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

% Compute \int_\Omega exp(k\pt(x-xc)) d \Omega
% xc is the center of the element


function f = integrate_exp_element(nodes,jk,xc)


if norm(jk)<1e-6
    f=polyarea(nodes(:,1),nodes(:,2));
else
    f=0;
    for i_node=2:size(nodes,1)
        [nx,ny]=normal_edge_out_element(nodes(i_node-1,:)',nodes(i_node,:)',xc);
        f=f-int_edge_1(jk,nodes(i_node-1,:)',nodes(i_node,:)',xc)*(jk(1)*nx+jk(2)*ny)/(norm(jk)^2);
    end
    [nx,ny]=normal_edge_out_element(nodes(i_node,:)',nodes(1,:)',xc);
    f=f-int_edge_1(jk,nodes(i_node,:)',nodes(1,:)',xc)*(jk(1)*nx+jk(2)*ny)/(norm(jk)^2)  ;
end


end



