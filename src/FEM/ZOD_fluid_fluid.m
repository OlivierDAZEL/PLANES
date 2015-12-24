% link_FEMZOD_fluid_fluid.m
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


for ii=1:size(edges.ZOD_moins,1)
    
    node_moins=edges.ZOD_moins(ii,[1 2 6])
    node_plus=edges.ZOD_plus(ii,[1 2 6])
    
    
    a1(1)=nodes(node_moins(1),1);
    a1(2)=nodes(node_moins(1),2);
    a2(1)=nodes(node_moins(2),1);
    a2(2)=nodes(node_moins(2),2);
    
    FSIe=TR6_FSI(a1,a2);
    
    
    index_force_p_moins=dof_A(p_TR(node_moins));
    index_F_elem_p_moins=find(index_force_p_moins);
    index_F_global_p_moins=index_force_p_moins(index_F_elem_p_moins);
    
    
    index_force_p_plus=dof_A(p_TR(node_plus));
    index_F_elem_p_plus=find(index_force_p_plus);
    index_F_global_p_plus=index_force_p_plus(index_F_elem_p_plus);
    
    A(index_F_global_p_moins,index_F_global_p_moins)=A(index_F_global_p_moins,index_F_global_p_moins)-TT(1,1)*(FSIe(index_F_elem_p_moins,index_F_elem_p_moins));
    A(index_F_global_p_moins,index_F_global_p_plus) =A(index_F_global_p_moins,index_F_global_p_plus) -TT(1,2)*(FSIe(index_F_elem_p_moins,index_F_elem_p_plus));
    A(index_F_global_p_plus,index_F_global_p_moins) =A(index_F_global_p_plus,index_F_global_p_moins) -TT(2,1)*(FSIe(index_F_elem_p_plus,index_F_elem_p_moins));
    A(index_F_global_p_plus,index_F_global_p_plus)  =A(index_F_global_p_plus,index_F_global_p_plus)  -TT(2,2)*(FSIe(index_F_elem_p_plus,index_F_elem_p_plus));
end