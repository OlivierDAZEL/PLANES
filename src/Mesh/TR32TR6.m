% TR32TR6.m
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


function [nb,nodes,elem,edges]=TR32TR6(nb_in,nodes_in,elem_in,edges_in)

% Creation of 3 segments by element
% segments=[node1 node2 #element1 0 element_label1]


nb=nb_in;
edges=edges_in;
elem=elem_in;

elements_TR6=find(elem_in.model==1);
nb_elements_TR6=length(elements_TR6);
elements_not_TR6=find(elem_in.model~=1);


node_supp=[(nodes_in(elem_in.nodes(elements_TR6,1),1)+nodes_in(elem_in.nodes(elements_TR6,2),1))/2 (nodes_in(elem_in.nodes(elements_TR6,1),2)+nodes_in(elem_in.nodes(elements_TR6,2),2))/2; ...
           (nodes_in(elem_in.nodes(elements_TR6,2),1)+nodes_in(elem_in.nodes(elements_TR6,3),1))/2 (nodes_in(elem_in.nodes(elements_TR6,2),2)+nodes_in(elem_in.nodes(elements_TR6,3),2))/2; ...
           (nodes_in(elem_in.nodes(elements_TR6,3),1)+nodes_in(elem_in.nodes(elements_TR6,1),1))/2 (nodes_in(elem_in.nodes(elements_TR6,3),2)+nodes_in(elem_in.nodes(elements_TR6,1),2))/2];

node_supp = unique(node_supp,'rows');

nodes=[nodes_in;node_supp];



elem.nodes(elements_not_TR6,1)=elem_in.nodes(elements_not_TR6,1);
elem.nodes(elements_not_TR6,2)=elem_in.nodes(elements_not_TR6,2);
elem.nodes(elements_not_TR6,3)=elem_in.nodes(elements_not_TR6,3);

elem.nodes(elements_TR6,1)=elem_in.nodes(elements_TR6,1);
elem.nodes(elements_TR6,3)=elem_in.nodes(elements_TR6,2);
elem.nodes(elements_TR6,5)=elem_in.nodes(elements_TR6,3);
nb.nodes=size(nodes,1);


for ie=1:nb_elements_TR6
    node_middle(nodes,elem.nodes(elements_TR6(ie),1),elem.nodes(elements_TR6(ie),3));
    elem.nodes(elements_TR6(ie),2)=node_middle(nodes,elem.nodes(elements_TR6(ie),1),elem.nodes(elements_TR6(ie),3));
    elem.nodes(elements_TR6(ie),4)=node_middle(nodes,elem.nodes(elements_TR6(ie),3),elem.nodes(elements_TR6(ie),5));
    elem.nodes(elements_TR6(ie),6)=node_middle(nodes,elem.nodes(elements_TR6(ie),5),elem.nodes(elements_TR6(ie),1));
end


if nb.dirichlets~=0
    for ie=1:nb.dirichlets
        if ismember(edges.dirichlets(ie,3),elements_TR6)
            edges.dirichlets(ie,6)=node_middle(nodes,edges.dirichlets(ie,1),edges.dirichlets(ie,2));
        end
    end
end
if nb.loads~=0
    for ie=1:nb.loads
        if ismember(edges.loads(ie,3),elements_TR6)
            edges.loads(ie,6)=node_middle(nodes,edges.loads(ie,1),edges.loads(ie,2));
        end
    end
end
if nb.periodicity~=0
    for ie=1:nb.periodicity
        if ismember(edges.periodicity(ie,3),elements_TR6)
            edges.periodicity(ie,6)=node_middle(nodes,edges.periodicity(ie,1),edges.periodicity(ie,2));
        end
    end
end
if nb.MMT~=0
    for ie=1:nb.MMT
        if ismember(edges.MMT(ie,3),elements_TR6)
            edges.MMT(ie,6)=node_middle(nodes,edges.MMT(ie,1),edges.MMT(ie,2));
        end
    end
end
if nb.internal~=0
    for ie=1:nb.internal
        if ((ismember(edges.internal(ie,3),elements_TR6))||(ismember(edges.internal(ie,4),elements_TR6)))
            edges.internal(ie,7)=node_middle(nodes,edges.internal(ie,1),edges.internal(ie,2));
        end
    end
end





end
