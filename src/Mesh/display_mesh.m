% display_mesh.m
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


% function void=display_mesh(nb,nodes,elements,element_label,edges,num_media,element_num_mat,interfaces,edges_MMT,loads,dirichlets,periodicity)

figure (10)
title('Nodes and Elements')
hold on
for ie=1:nb.elements
    if elem.model(ie)==1
        line([nodes(elem.nodes(ie,1),1) nodes(elem.nodes(ie,3),1)],[nodes(elem.nodes(ie,1),2) nodes(elem.nodes(ie,3),2)],'Color','r');
        line([nodes(elem.nodes(ie,3),1) nodes(elem.nodes(ie,5),1)],[nodes(elem.nodes(ie,3),2) nodes(elem.nodes(ie,5),2)],'Color','r');
        line([nodes(elem.nodes(ie,5),1) nodes(elem.nodes(ie,1),1)],[nodes(elem.nodes(ie,5),2) nodes(elem.nodes(ie,1),2)],'Color','r');
        text(mean(nodes(elem.nodes(ie,1:6),1)),mean(nodes(elem.nodes(ie,1:6),2)),num2str(ie),'Fontsize',15);
    end
    if elem.model(ie)==2
        line([nodes(elem.nodes(ie,1),1) nodes(elem.nodes(ie,2),1)],[nodes(elem.nodes(ie,1),2) nodes(elem.nodes(ie,2),2)],'Color','r');
        line([nodes(elem.nodes(ie,2),1) nodes(elem.nodes(ie,3),1)],[nodes(elem.nodes(ie,2),2) nodes(elem.nodes(ie,3),2)],'Color','r');
        line([nodes(elem.nodes(ie,3),1) nodes(elem.nodes(ie,4),1)],[nodes(elem.nodes(ie,3),2) nodes(elem.nodes(ie,4),2)],'Color','r');
        line([nodes(elem.nodes(ie,4),1) nodes(elem.nodes(ie,1),1)],[nodes(elem.nodes(ie,4),2) nodes(elem.nodes(ie,1),2)],'Color','r');
        text(mean(nodes(elem.nodes(ie,1:4),1)),mean(nodes(elem.nodes(ie,1:4),2)),num2str(ie),'Fontsize',15);
    end
    
    if ismember(elem.model(ie),[3,10])
        line([nodes(elem.nodes(ie,1),1) nodes(elem.nodes(ie,2),1)],[nodes(elem.nodes(ie,1),2) nodes(elem.nodes(ie,2),2)],'Color','r');
        line([nodes(elem.nodes(ie,2),1) nodes(elem.nodes(ie,3),1)],[nodes(elem.nodes(ie,2),2) nodes(elem.nodes(ie,3),2)],'Color','r');
        line([nodes(elem.nodes(ie,3),1) nodes(elem.nodes(ie,1),1)],[nodes(elem.nodes(ie,3),2) nodes(elem.nodes(ie,1),2)],'Color','r');
        text(mean(nodes(elem.nodes(ie,1:3),1)),mean(nodes(elem.nodes(ie,1:3),2)),num2str(ie),'Fontsize',15);
    end
end
plot(nodes(:,1),nodes(:,2),'b.','Markersize',15);
for ii=1:nb.nodes
    text(nodes(ii,1),nodes(ii,2),num2str(ii),'Fontsize',15);
end
axis equal




figure(11)
title('Edges and Models')
hold on
for ie=1:nb.elements
    if elem.model(ie)==1
        line([nodes(elem.nodes(ie,1),1) nodes(elem.nodes(ie,3),1)],[nodes(elem.nodes(ie,1),2) nodes(elem.nodes(ie,3),2)],'Color','r');
        line([nodes(elem.nodes(ie,3),1) nodes(elem.nodes(ie,5),1)],[nodes(elem.nodes(ie,3),2) nodes(elem.nodes(ie,5),2)],'Color','r');
        line([nodes(elem.nodes(ie,5),1) nodes(elem.nodes(ie,1),1)],[nodes(elem.nodes(ie,5),2) nodes(elem.nodes(ie,1),2)],'Color','r');
        text(mean(nodes(elem.nodes(ie,1:6),1)),mean(nodes(elem.nodes(ie,1:6),2)),num2str(elem.model(ie)),'Fontsize',15);
    end
    if elem.model(ie)==2
        line([nodes(elem.nodes(ie,1),1) nodes(elem.nodes(ie,2),1)],[nodes(elem.nodes(ie,1),2) nodes(elem.nodes(ie,2),2)],'Color','r');
        line([nodes(elem.nodes(ie,2),1) nodes(elem.nodes(ie,3),1)],[nodes(elem.nodes(ie,2),2) nodes(elem.nodes(ie,3),2)],'Color','r');
        line([nodes(elem.nodes(ie,3),1) nodes(elem.nodes(ie,4),1)],[nodes(elem.nodes(ie,3),2) nodes(elem.nodes(ie,4),2)],'Color','r');
        line([nodes(elem.nodes(ie,4),1) nodes(elem.nodes(ie,1),1)],[nodes(elem.nodes(ie,4),2) nodes(elem.nodes(ie,1),2)],'Color','r');
        text(mean(nodes(elem.nodes(ie,1:4),1)),mean(nodes(elem.nodes(ie,1:4),2)),num2str(elem.model(ie)),'Fontsize',15);
    end
    if ismember(elem.model(ie),[3,10])
        line([nodes(elem.nodes(ie,1),1) nodes(elem.nodes(ie,2),1)],[nodes(elem.nodes(ie,1),2) nodes(elem.nodes(ie,2),2)],'Color','r');
        line([nodes(elem.nodes(ie,2),1) nodes(elem.nodes(ie,3),1)],[nodes(elem.nodes(ie,2),2) nodes(elem.nodes(ie,3),2)],'Color','r');
        line([nodes(elem.nodes(ie,3),1) nodes(elem.nodes(ie,1),1)],[nodes(elem.nodes(ie,3),2) nodes(elem.nodes(ie,1),2)],'Color','r');
        text(mean(nodes(elem.nodes(ie,1:3),1)),mean(nodes(elem.nodes(ie,1:3),2)),num2str(elem.model(ie)),'Fontsize',15);
    end
end

for ii=1:nb.loads
    line([nodes(edges.loads(ii,1),1) nodes(edges.loads(ii,2),1)],[nodes(edges.loads(ii,1),2) nodes(edges.loads(ii,2),2)]);
    text((nodes(edges.loads(ii,1),1)+nodes(edges.loads(ii,2),1))/2,(nodes(edges.loads(ii,1),2)+nodes(edges.loads(ii,2),2))/2,num2str(edges.loads(ii,4)),'Fontsize',15);
end
for ii=1:nb.periodicity
    line([nodes(edges.periodicity(ii,1),1) nodes(edges.periodicity(ii,2),1)],[nodes(edges.periodicity(ii,1),2) nodes(edges.periodicity(ii,2),2)]);
    text((nodes(edges.periodicity(ii,1),1)+nodes(edges.periodicity(ii,2),1))/2,(nodes(edges.periodicity(ii,1),2)+nodes(edges.periodicity(ii,2),2))/2,num2str(edges.periodicity(ii,4)),'Fontsize',15);
    
end
for ii=1:nb.MMT
    line([nodes(edges_MMT(ii,1),1) nodes(edges_MMT(ii,2),1)],[nodes(edges_MMT(ii,1),2) nodes(edges_MMT(ii,2),2)]);
    text((nodes(edges_MMT(ii,1),1)+nodes(edges_MMT(ii,2),1))/2,(nodes(edges_MMT(ii,1),2)+nodes(edges_MMT(ii,2),2))/2,num2str(edges_MMT(ii,4)),'Fontsize',15);
    
end
for ii=1:nb.dirichlets
    line([nodes(edges.dirichlets(ii,1),1) nodes(edges.dirichlets(ii,2),1)],[nodes(edges.dirichlets(ii,1),2) nodes(edges.dirichlets(ii,2),2)]);
    text((nodes(edges.dirichlets(ii,1),1)+nodes(edges.dirichlets(ii,2),1))/2,(nodes(edges.dirichlets(ii,1),2)+nodes(edges.dirichlets(ii,2),2))/2,num2str(edges.dirichlets(ii,4)),'Fontsize',15);
    
end
for ii=1:nb.internal
    line([nodes(edges.internal(ii,1),1) nodes(edges.internal(ii,2),1)],[nodes(edges.internal(ii,1),2) nodes(edges.internal(ii,2),2)]);
    text((nodes(edges.internal(ii,1),1)+nodes(edges.internal(ii,2),1))/2,(nodes(edges.internal(ii,1),2)+nodes(edges.internal(ii,2),2))/2,num2str(edges.internal(ii,4)),'Fontsize',15);
    
end
axis equal


figure (13)
title('Materials')
hold on
for ie=1:nb.elements
    if elem.model(ie)==1
        line([nodes(elem.nodes(ie,1),1) nodes(elem.nodes(ie,3),1)],[nodes(elem.nodes(ie,1),2) nodes(elem.nodes(ie,3),2)],'Color','r');
        line([nodes(elem.nodes(ie,3),1) nodes(elem.nodes(ie,5),1)],[nodes(elem.nodes(ie,3),2) nodes(elem.nodes(ie,5),2)],'Color','r');
        line([nodes(elem.nodes(ie,5),1) nodes(elem.nodes(ie,1),1)],[nodes(elem.nodes(ie,5),2) nodes(elem.nodes(ie,1),2)],'Color','r');
        text(mean(nodes(elem.nodes(ie,1:6),1)),mean(nodes(elem.nodes(ie,1:6),2)),num2str(elem.label(ie)),'Fontsize',15);
    end
        if elem.model(ie)==2
        line([nodes(elem.nodes(ie,1),1) nodes(elem.nodes(ie,2),1)],[nodes(elem.nodes(ie,1),2) nodes(elem.nodes(ie,2),2)],'Color','r');
        line([nodes(elem.nodes(ie,2),1) nodes(elem.nodes(ie,3),1)],[nodes(elem.nodes(ie,2),2) nodes(elem.nodes(ie,3),2)],'Color','r');
        line([nodes(elem.nodes(ie,3),1) nodes(elem.nodes(ie,4),1)],[nodes(elem.nodes(ie,3),2) nodes(elem.nodes(ie,4),2)],'Color','r');
        line([nodes(elem.nodes(ie,4),1) nodes(elem.nodes(ie,1),1)],[nodes(elem.nodes(ie,4),2) nodes(elem.nodes(ie,1),2)],'Color','r');
        text(mean(nodes(elem.nodes(ie,1:4),1)),mean(nodes(elem.nodes(ie,1:4),2)),num2str(elem.label(ie)),'Fontsize',15);
    end
    if ismember(elem.model(ie),[3,10])
        line([nodes(elem.nodes(ie,1),1) nodes(elem.nodes(ie,2),1)],[nodes(elem.nodes(ie,1),2) nodes(elem.nodes(ie,2),2)],'Color','r');
        line([nodes(elem.nodes(ie,2),1) nodes(elem.nodes(ie,3),1)],[nodes(elem.nodes(ie,2),2) nodes(elem.nodes(ie,3),2)],'Color','r');
        line([nodes(elem.nodes(ie,3),1) nodes(elem.nodes(ie,1),1)],[nodes(elem.nodes(ie,3),2) nodes(elem.nodes(ie,1),2)],'Color','r');
        text(mean(nodes(elem.nodes(ie,1:3),1)),mean(nodes(elem.nodes(ie,1:3),2)),num2str(elem.label(ie)),'Fontsize',15);
    end
end
axis equal

% end
