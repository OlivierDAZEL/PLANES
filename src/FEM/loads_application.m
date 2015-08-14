% loads_application_TR6.m
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

for ie=1:nb.loads
    typ=floor(edges.loads(ie,4));
    
    x1=nodes(edges.loads(ie,1),1);
    y1=nodes(edges.loads(ie,1),2);
    x2=nodes(edges.loads(ie,2),1);
    y2=nodes(edges.loads(ie,2),2);
    length_edge=sqrt((x2-x1)^2+(y2-y1)^2);
    
    switch typ
        case {3}  % Unit
            
            if elem.model(edges.loads(ie,4))==1
                if (x1<x2)
                    a=x1;
                    node(1)=edges.loads(ie,1);
                    node(2)=edges.loads(ie,2);
                    node(3)=edges.loads(ie,6);
                else
                    a=x2;
                    node(2)=edges.loads(ie,1);
                    node(1)=edges.loads(ie,2);
                    node(3)=edges.loads(ie,6);
                end
                F3=TR6_unit(length_edge);
                index_force=dof_A(p_TR(node));
                index_F_elem=find(index_force);
                index_F_global=index_force(index_F_elem);
                F(index_F_global)=F(index_F_global)-F3(index_F_elem)/(1j*omega);
            elseif  elem.model(edges.loads(ie,4))==2
                lx=norm(nodes(elem.nodes(edges.loads(ie,3),1),:)-nodes(elem.nodes(edges.loads(ie,3),2),:));
                ly=norm(nodes(elem.nodes(edges.loads(ie,3),1),:)-nodes(elem.nodes(edges.loads(ie,3),4),:));
                load_Hermite_2D_2
                index_p=dof_A(p_H12(elem.nodes(edges.loads(ie,3),1:4)));
                % What is the edge on the element ?
                if sort(loads(ie,1:2))==sort(elements(loads(ie,3),1:2)) % Bottom
                    indice_test=index_p([1 2 5 4]);
                elseif sort(loads(ie,1:2))==sort(elements(loads(ie,3),2:3))
                    indice_test=index_p([4 6 9 7]);
                elseif sort(loads(ie,1:2))==sort(elements(loads(ie,3),2:3))
                    indice_test=index_p([10 11 8 7]);
                elseif sort(loads(ie,1:2))==sort(elements(loads(ie,3),2:3))
                    indice_test=index_p([1 3 12 10]);
                end
                for i_test=1:4
                    eval(['Psi_test=Psi_',num2str(i_test),'_x;'])
                    F(indice_test(i_test))=F(indice_test(i_test))+integrate_polynom(Psi_test,length_edge)/(1j*omega);
                end
            end
    end
end
