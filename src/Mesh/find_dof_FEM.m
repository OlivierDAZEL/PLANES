% find_dof_FEM.m
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


isvalidof_FEM=zeros(12*nb.nodes,1);
dof_A=zeros(12*nb.nodes,1);



for ie=1:nb.elements
    typ=floor(elem.label(ie)/1000);
    if ismember(elem.model(ie),[1 3])
        switch typ
            case {0,2,3,8} %! Acoustic/EF/limp/PML
                isvalidof_FEM(p_TR(nonzeros(elem.nodes(ie,1:6))))=1;
            case {1} %! Elastic solid
                isvalidof_FEM(uxy_TR(nonzeros(elem.nodes(ie,1:6))))=1;
            case {4,5}	%! PEM
                isvalidof_FEM(uxyp_TR(nonzeros(elem.nodes(ie,1:6))))=1;
            otherwise
                disp('Attention element sans nature connue')
                stop
        end
    end
        if ismember(elem.model(ie),2)
        switch typ
            case {0,2,3,8} %! Acoustic/EF/limp/PML
                isvalidof_FEM(p_H12(nonzeros(elem.nodes(ie,1:4))))=1;
            case {1} %! Elastic solid
                isvalidof_FEM(uxy_H12(nonzeros(elem.nodes(ie,1:4))))=1;
            case {4,5}	%! PEM
                isvalidof_FEM(uxyp_12(nonzeros(elem.nodes(ie,1:4))))=1;
            otherwise
                disp('Attention element sans nature connue')
                stop
        end
    end
 end



for ie=1:nb.dirichlets
    if (edges.dirichlets(ie,4)==5) % Sliding
        xx=abs(nodes(edges.dirichlets(ie,1),1)-nodes(edges.dirichlets(ie,2),1));
        yy=abs(nodes(edges.dirichlets(ie,1),2)-nodes(edges.dirichlets(ie,2),2));
        if (xx>yy)
            isvalidof_FEM(uy_H16(edges.dirichlets(ie,:)))=0;
        else
            isvalidof_FEM(ux_H16(edges.dirichlets(ie,:)))=0;
        end
    end
    if (edges.dirichlets(ie,4)==6) % Bonded
            isvalidof_FEM(uxy_H16(edges.dirichlets(ie,:)))=0;
    end
end



dof_back=[];

itemp=0;
for ii=1:12*nb.nodes
    if (isvalidof_FEM(ii))
        itemp=itemp+1;
        dof_A(ii)=itemp;
        dof_back(itemp)=ii;
    end
end

nb.dof_FEM=itemp;
list_dof_valid=find(isvalidof_FEM);
clear isvalidof_FEM

