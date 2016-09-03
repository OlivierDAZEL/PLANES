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


function void=display_mesh_light(nb,nodes,elem,Linewidth)

hold on
for ie=1:nb.elements
    if elem.model(ie)==1
        line([nodes(elem.nodes(ie,1),1) nodes(elem.nodes(ie,3),1)],[nodes(elem.nodes(ie,1),2) nodes(elem.nodes(ie,3),2)],'Color','k','Linewidth',Linewidth);
        line([nodes(elem.nodes(ie,3),1) nodes(elem.nodes(ie,5),1)],[nodes(elem.nodes(ie,3),2) nodes(elem.nodes(ie,5),2)],'Color','k','Linewidth',Linewidth);
        line([nodes(elem.nodes(ie,5),1) nodes(elem.nodes(ie,1),1)],[nodes(elem.nodes(ie,5),2) nodes(elem.nodes(ie,1),2)],'Color','k','Linewidth',Linewidth);
    end
    if ismember(elem.model(ie),[2,11])
        line([nodes(elem.nodes(ie,1),1) nodes(elem.nodes(ie,2),1)],[nodes(elem.nodes(ie,1),2) nodes(elem.nodes(ie,2),2)],'Color','k','Linewidth',Linewidth);
        line([nodes(elem.nodes(ie,2),1) nodes(elem.nodes(ie,3),1)],[nodes(elem.nodes(ie,2),2) nodes(elem.nodes(ie,3),2)],'Color','k','Linewidth',Linewidth);
        line([nodes(elem.nodes(ie,3),1) nodes(elem.nodes(ie,4),1)],[nodes(elem.nodes(ie,3),2) nodes(elem.nodes(ie,4),2)],'Color','k','Linewidth',Linewidth);
        line([nodes(elem.nodes(ie,4),1) nodes(elem.nodes(ie,1),1)],[nodes(elem.nodes(ie,4),2) nodes(elem.nodes(ie,1),2)],'Color','k','Linewidth',Linewidth);
    end
    
    if ismember(elem.model(ie),[3,10])
        line([nodes(elem.nodes(ie,1),1) nodes(elem.nodes(ie,2),1)],[nodes(elem.nodes(ie,1),2) nodes(elem.nodes(ie,2),2)],'Color','k','Linewidth',Linewidth);
        line([nodes(elem.nodes(ie,2),1) nodes(elem.nodes(ie,3),1)],[nodes(elem.nodes(ie,2),2) nodes(elem.nodes(ie,3),2)],'Color','k','Linewidth',Linewidth);
        line([nodes(elem.nodes(ie,3),1) nodes(elem.nodes(ie,1),1)],[nodes(elem.nodes(ie,3),2) nodes(elem.nodes(ie,1),2)],'Color','k','Linewidth',Linewidth);
    end
end

