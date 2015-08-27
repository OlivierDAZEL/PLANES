% createmshH12.m
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


function [nb,nodes,elem,edges]=createmshH12(lx,ly,nx,ny,label_boundaries)

nb.nodes=(nx+1)*(ny+1);
nb.elements=nx*ny;
nb.edges=2*(nx+ny);

nodes=zeros(nb.nodes,2);
elem.nodes=zeros(nb.elements,3);
elem.label=zeros(nb.elements,1);
elem.model=zeros(nb.elements,1);
edges=zeros(nb.edges,3);

temp=0;
Delta_x=lx/nx;
Delta_y=ly/ny;
for iy=1:ny+1
    ynode=(iy-1)*Delta_y;
    
    for ix=1:nx+1
        xnode=(ix-1)*Delta_x;
        temp=temp+1;
        nodes(temp,1)=xnode;
        nodes(temp,2)=ynode;
    end
end

temp=0;
for ix=1:nx
    for iy=1:ny
        temp=temp+1;
        node_base=ix+(iy-1)*(nx+1);
        elem.nodes(temp,1)=node_base;
        elem.nodes(temp,2)=node_base+1;
        elem.nodes(temp,3)=node_base+(nx+2);
        elem.nodes(temp,4)=node_base+(nx+1);
        elem.label(temp)=0;
    end
end




temp=0;
for ix=1:nx
    temp=temp+1;
    node_base=ix;
    edges(temp,1)=node_base;
    edges(temp,2)=node_base+1;
    edges(temp,3)=label_boundaries(1);
end

for ix=1:nx
    temp=temp+1;
    node_base=ix+ny*(nx+1);
    edges(temp,1)=node_base;
    edges(temp,2)=node_base+1;
    edges(temp,3)=label_boundaries(3);
end

for iy=1:ny
    temp=temp+1;
    node_base=nx+1+(iy-1)*(nx+1);
    edges(temp,1)=node_base;
    edges(temp,2)=node_base+nx+1;
    edges(temp,3)=label_boundaries(2);
end

for iy=1:ny
    temp=temp+1;
    node_base=1+(iy-1)*(nx+1);
    edges(temp,1)=node_base;
    edges(temp,2)=node_base+nx+1;
    edges(temp,3)=label_boundaries(4);
end


