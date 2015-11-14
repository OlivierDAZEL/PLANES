% msh_import.m
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


function [nb,nodes,elem,edge_msh]=msh_import(nom_fichier)

fid=fopen(nom_fichier,'r');
nb.nodes=fscanf(fid,'%i',1);
nb.elements=fscanf(fid,'%i',1);
nb.edges=fscanf(fid,'%i',1);

nodes=zeros(nb.nodes,2);
elem.nodes=zeros(nb.elements,3);
elem.label=zeros(nb.elements,1);
edge_msh=zeros(nb.edges,3);


for ii=1:nb.nodes
    nodes(ii,1)=fscanf(fid,'%f',1);
    nodes(ii,2)=fscanf(fid,'%f',1);
    node_label(ii)=fscanf(fid,'%f',1);
end
for ii=1:nb.elements
    elem.nodes(ii,1)=fscanf(fid,'%i',1);
    elem.nodes(ii,2)=fscanf(fid,'%i',1);
    elem.nodes(ii,3)=fscanf(fid,'%i',1);
    elem.label(ii,1)=abs(fscanf(fid,'%f',1));
    
end
for ii=1:nb.edges
    edge_msh(ii,1)=fscanf(fid,'%i',1);
    edge_msh(ii,2)=fscanf(fid,'%i',1);
    edge_msh(ii,3)=fscanf(fid,'%i',1);
end

fclose(fid);


