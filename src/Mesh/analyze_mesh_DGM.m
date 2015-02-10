% analyze_mesh_DGM.m
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


period=max(nodes(:,1))-min(nodes(:,1));

ondes_element=zeros(nb_elements,1);
for ie=1:nb_elements
   switch element_label(ie)
       case {0,2000:2999,3000:3999}
           ondes_element(ie)=1;
       case {1000:1999}
           ondes_element(ie)=2;   
       case {4000:5999}
           ondes_element(ie)=3;
   end
end
dof_start_element=zeros(nb_elements,1);
dof_start_element(1)=1;
for ie=2:nb_elements
   dof_start_element(ie)=dof_start_element(ie-1)+ondes_element(ie-1)*nb_theta; 
end


nb_dof_DGM=dof_start_element(ie)+ondes_element(ie)*nb_theta-1; 

vec_theta=linspace(0,2*pi,nb_theta+1);
vec_theta(end)=[];


if nb_periodicity~=0
    
    edge_left= find(periodicity(:,4)==98);
    edge_right=find(periodicity(:,4)==99);
    
    y_left=sort([nodes(periodicity(edge_left,1),2) nodes(periodicity(edge_left,2),2)],2);
    y_right=sort([nodes(periodicity(edge_right,1),2) nodes(periodicity(edge_right,2),2)],2);

    qdsqsddsq
    
    [temp,i_left]=sort(nodes(node_left,2));
    node_left=node_left(i_left);
    
    node_right=unique([periodicity(edge_right,1);periodicity(edge_right,2)]);
    [temp,i_right]=sort(nodes(node_right,2));
    node_right=node_right(i_right);
    
    
end


rzeerzezrez




