% apply_FSI.m
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


for ie=1:nb.internal

    x1=nodes(edges.internal(ie,1),1);
    y1=nodes(edges.internal(ie,1),2);
    x2=nodes(edges.internal(ie,2),1);
    y2=nodes(edges.internal(ie,2),2);
    length_edge=sqrt((x2-x1)^2+(y2-y1)^2);
    if (x1<x2)
        a=x1;
        node(1)=edges.internal(ie,1);
        node(2)=edges.internal(ie,2);
        node(3)=edges.internal(ie,7);
    else
        a=x2;
        node(2)=edges.internal(ie,1);
        node(1)=edges.internal(ie,2);
        node(3)=edges.internal(ie,7);
    end
    
    
    a1(1)=nodes(node(1),1);
    a1(2)=nodes(node(1),2);
    a2(1)=nodes(node(2),1);
    a2(2)=nodes(node(2),2);
    
    vec_tangent=a2-a1;
    vec_normal=[vec_tangent(2) -vec_tangent(1)];
    vec_normal=vec_normal/norm(vec_normal); 
    num_element=edges.internal(ie,3);
    center_e1=mean(nodes(nonzeros(elem.nodes(num_element,:)),1:2))';

    
    vec_temp=nodes(edges.internal(ie,7),1:2)-center_e1';
    temp=vec_normal*vec_temp';
    if temp<0
        vec_normal=-vec_normal;
    end
     
    if floor(elem.label(edges.internal(ie,3))/1000)==0
        n_elas=-vec_normal;
        n_acou= n_elas;
    else
        n_elas= vec_normal;
        n_acou= n_elas;
    end
    n_elas=-n_elas;
    n_acou=-n_acou;

    
    FSIe=TR6_FSI(a1,a2);
    
    index_force_p=dof_A(p_TR(node));
    index_F_elem_p=find(index_force_p);
    index_F_global_p=index_force_p(index_F_elem_p);
    
    
    index_force_ux=dof_A(ux_TR(node));
    index_F_elem_ux=find(index_force_ux);
    index_F_global_ux=index_force_ux(index_F_elem_ux);
   
    index_force_uy=dof_A(uy_TR(node));
    index_F_elem_uy=find(index_force_uy);
    index_F_global_uy=index_force_uy(index_F_elem_uy);


    A(index_F_global_p,index_F_global_ux)=A(index_F_global_p,index_F_global_ux)-n_acou(1)*FSIe(index_F_elem_p,index_F_elem_ux);    
    A(index_F_global_p,index_F_global_uy)=A(index_F_global_p,index_F_global_uy)-n_acou(2)*FSIe(index_F_elem_p,index_F_elem_uy);
    A(index_F_global_ux,index_F_global_p)=A(index_F_global_ux,index_F_global_p)-n_elas(1)*FSIe(index_F_elem_ux,index_F_elem_p);    
    A(index_F_global_uy,index_F_global_p)=A(index_F_global_uy,index_F_global_p)-n_elas(2)*FSIe(index_F_elem_uy,index_F_elem_p);    
    
    
end