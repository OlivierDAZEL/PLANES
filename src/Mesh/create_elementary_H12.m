% create_elementary_H12.m
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


function [elem,H_elem_H12,Q_elem_H12]=create_elementary_H12(nb_in,nodes_in,elem_in)

% Find the number of H12


elem=elem_in;

index_H12=find(elem_in.model==2);

for ii=1:length(index_H12)
    lx(ii)=norm(nodes_in(elem_in.nodes(index_H12(ii),1),:)-nodes_in(elem_in.nodes(index_H12(ii),2),:));
    ly(ii)=norm(nodes_in(elem_in.nodes(index_H12(ii),1),:)-nodes_in(elem_in.nodes(index_H12(ii),4),:));
end


% For roundoff errors
lc=1e8*(lx+1j*ly);
lc=round(lc)/1e8;

[lc,~,index] = unique(lc);

H_elem_H12=zeros(12,12,length(lc));
Q_elem_H12=zeros(12,12,length(lc));

for ii=1:length(lc)
    [H_elem_H12(1:12,1:12,ii),Q_elem_H12(1:12,1:12,ii)] = H12_fluid(real(lc(ii)),imag(lc(ii)));
end
elem.H12(index_H12)=index;
end
