% edge_periodicity.m
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


%%%%% Coordinates of the edge's vertices on the right and left edges

coord_edge_left(1:2,1)=nodes(periodicity_left(ie,1),1:2)';
coord_edge_left(1:2,2)=nodes(periodicity_left(ie,2),1:2)';

a_left=coord_edge_left(:,1);
b_left=coord_edge_left(:,2);

coord_edge_right(1:2,1)=nodes(periodicity_right(ie,1),1:2)';
coord_edge_right(1:2,2)=nodes(periodicity_right(ie,2),1:2)';

a_right=coord_edge_right(:,1);
b_right=coord_edge_right(:,2);




%%%%% Elements on both sides of the edge

e_right=periodicity_right(ie,3);
e_left =periodicity_left(ie,3);

c_left=mean(nodes(elements(e_left,:),1:2))';
c_right=mean(nodes(elements(e_right,:),1:2))';


%%%%% vector normal pointing outwards \Omega_right

nx=-1;
ny=0;

e_edge=e_right;
parameter_element


valid_edge_internal=0;
if (element_label(e_1)==element_label(e_2))
    if ((floor(element_label(e_edge)/1000)==0)|(floor(element_label(e_edge)/1000)==2)|(floor(element_label(e_edge)/1000)==3))
        %disp('Lancement internal_fluid')
        fluid_periodic
        valid_edge_internal=1;
     elseif (floor(element_label(e_edge)/1000)==4)
%         %disp('Lancement internal_PEM')
         PEM_periodic
         valid_edge_internal=1;
%     elseif (floor(element_label(e_edge)/1000)==8)
%         parameter_element
%         internal_PML
%         valid_edge_internal=1;
     end
%     
    
else % (element_label(e_1)~=element_label(e_2))
    
    
%     
%     
%     if (sum(floor(element_label(e_1)/1000)==[0 2 3]))*(sum(floor(element_label(e_2)/1000)==[0 2 3]))
%         %disp('Lancement fluid_fluid')
%         
%         fluid_fluid
%         valid_edge_internal=1;
%     end
%     if (sum(floor(element_label(e_1)/1000)==[0 2 3]))*(sum(floor(element_label(e_2)/1000)==[4]))
%         %disp('Lancement fluid_PEM')
%         fluid_PEM
%         valid_edge_internal=1;
%     end
%     if (sum(floor(element_label(e_2)/1000)==[0 2 3]))*(sum(floor(element_label(e_1)/1000)==[4]))
%         %disp('Lancement PEM_fluid')
%         PEM_fluid
%         valid_edge_internal=1;
%     end
%     
%         if (sum(floor(element_label(e_1)/1000)==[0 2 3 8]))*(sum(floor(element_label(e_2)/1000)==[8]))
%         %disp('Lancement PML_PML')
%         PML_PML
%         valid_edge_internal=1;
%     end

end

if valid_edge_internal==0
    disp('Stop in edge internal not a valid internal edge')
    jkljkkjklj
end


