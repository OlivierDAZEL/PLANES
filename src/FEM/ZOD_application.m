% link_FEMZOD_elas_elas.m
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
    
    number_ZOD=1+(edges.ZOD(index_ZOD_moins(ii),4)-401)/2;
    TT=build_FEM_transfer_2D(k_air*sin(data_model.multilayer_ZOD(number_ZOD).theta_ZOD),elem.label(edges.ZOD(ii,3)),elem.label(edges.ZOD(ii,3)),omega,data_model.multilayer_ZOD(number_ZOD),k_air,air);

    
    node_moins=edges.ZOD_moins(ii,[1 2 6]);
    node_plus=edges.ZOD_plus(ii,[1 2 6]);
    
    a1(1)=nodes(node_moins(1),1);
    a1(2)=nodes(node_moins(1),2);
    a2(1)=nodes(node_moins(2),1);
    a2(2)=nodes(node_moins(2),2);
    
    switch floor(elem.label(edges.ZOD_moins(ii,3))/1000)
        
        case {0}
            
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
            
            
            
        case {1}
                        
            
            vec_tangent=a2-a1;
            vec_normal=[vec_tangent(2) vec_tangent(1)];
            vec_normal=vec_normal/norm(vec_normal);
            
            center_element=mean(nodes(elem.nodes(edges.ZOD_moins(ii,3),:),1:2),1);
            vec_temp=nodes(node_moins(3),1:2)-center_element;
            temp=vec_normal*vec_temp';
            if temp<0
                vec_normal=-vec_normal;
            end
            theta_rot=angle(vec_normal(1)+1i*vec_normal(2))-pi/2;
            Mat_rot=[cos(theta_rot) sin(theta_rot);-sin(theta_rot) cos(theta_rot)];
            Mat_rot=[Mat_rot 0*Mat_rot;0*Mat_rot Mat_rot];
            Mat_rotm1=[cos(theta_rot) -sin(theta_rot);sin(theta_rot) cos(theta_rot)];
            Mat_rotm1=[Mat_rotm1 0*Mat_rotm1;0*Mat_rotm1 Mat_rotm1];
            TTrot=Mat_rotm1*TT*Mat_rot;
            

            
            FSIe=TR6_FSI(a1,a2);
            
            index_force_ux_moins=dof_A(ux_TR(node_moins));
            index_F_elem_ux_moins=find(index_force_ux_moins);
            index_F_global_ux_moins=index_force_ux_moins(index_F_elem_ux_moins);
            
            index_force_uy_moins=dof_A(uy_TR(node_moins));
            index_F_elem_uy_moins=find(index_force_uy_moins);
            index_F_global_uy_moins=index_force_uy_moins(index_F_elem_uy_moins);
            
            index_force_ux_plus=dof_A(ux_TR(node_plus));
            index_F_elem_ux_plus=find(index_force_ux_plus);
            index_F_global_ux_plus=index_force_ux_plus(index_F_elem_ux_plus);
            
            
            index_force_uy_plus=dof_A(uy_TR(node_plus));
            index_F_elem_uy_plus=find(index_force_uy_plus);
            index_F_global_uy_plus=index_force_uy_plus(index_F_elem_uy_plus);
            
            
            A(index_F_global_ux_moins,index_F_global_ux_moins)=A(index_F_global_ux_moins,index_F_global_ux_moins)-TTrot(1,1)*(FSIe(index_F_elem_ux_moins,index_F_elem_ux_moins));
            A(index_F_global_ux_moins,index_F_global_uy_moins)=A(index_F_global_ux_moins,index_F_global_uy_moins)-TTrot(1,2)*(FSIe(index_F_elem_ux_moins,index_F_elem_uy_moins));
            A(index_F_global_ux_moins,index_F_global_ux_plus )=A(index_F_global_ux_moins,index_F_global_ux_plus )-TTrot(1,3)*(FSIe(index_F_elem_ux_moins,index_F_elem_ux_plus ));
            A(index_F_global_ux_moins,index_F_global_uy_plus )=A(index_F_global_ux_moins,index_F_global_uy_plus )-TTrot(1,4)*(FSIe(index_F_elem_ux_moins,index_F_elem_uy_plus ));
            
            A(index_F_global_uy_moins,index_F_global_ux_moins)=A(index_F_global_uy_moins,index_F_global_ux_moins)-TTrot(2,1)*(FSIe(index_F_elem_uy_moins,index_F_elem_ux_moins));
            A(index_F_global_uy_moins,index_F_global_uy_moins)=A(index_F_global_uy_moins,index_F_global_uy_moins)-TTrot(2,2)*(FSIe(index_F_elem_uy_moins,index_F_elem_uy_moins));
            A(index_F_global_uy_moins,index_F_global_ux_plus )=A(index_F_global_uy_moins,index_F_global_ux_plus )-TTrot(2,3)*(FSIe(index_F_elem_uy_moins,index_F_elem_ux_plus ));
            A(index_F_global_uy_moins,index_F_global_uy_plus )=A(index_F_global_uy_moins,index_F_global_uy_plus )-TTrot(2,4)*(FSIe(index_F_elem_uy_moins,index_F_elem_uy_plus ));
            
            
            A(index_F_global_ux_plus ,index_F_global_ux_moins)=A(index_F_global_ux_plus ,index_F_global_ux_moins)-TTrot(3,1)*(FSIe(index_F_elem_ux_plus ,index_F_elem_ux_moins));
            A(index_F_global_ux_plus ,index_F_global_uy_moins)=A(index_F_global_ux_plus ,index_F_global_uy_moins)-TTrot(3,2)*(FSIe(index_F_elem_ux_plus ,index_F_elem_uy_moins));
            A(index_F_global_ux_plus ,index_F_global_ux_plus )=A(index_F_global_ux_plus ,index_F_global_ux_plus )-TTrot(3,3)*(FSIe(index_F_elem_ux_plus ,index_F_elem_ux_plus ));
            A(index_F_global_ux_plus ,index_F_global_uy_plus )=A(index_F_global_ux_plus ,index_F_global_uy_plus )-TTrot(3,4)*(FSIe(index_F_elem_ux_plus ,index_F_elem_uy_plus ));
            
            A(index_F_global_uy_plus ,index_F_global_ux_moins)=A(index_F_global_uy_plus ,index_F_global_ux_moins)-TTrot(4,1)*(FSIe(index_F_elem_uy_plus ,index_F_elem_ux_moins));
            A(index_F_global_uy_plus ,index_F_global_uy_moins)=A(index_F_global_uy_plus ,index_F_global_uy_moins)-TTrot(4,2)*(FSIe(index_F_elem_uy_plus ,index_F_elem_uy_moins));
            A(index_F_global_uy_plus ,index_F_global_ux_plus )=A(index_F_global_uy_plus ,index_F_global_ux_plus )-TTrot(4,3)*(FSIe(index_F_elem_uy_plus ,index_F_elem_ux_plus ));
            A(index_F_global_uy_plus ,index_F_global_uy_plus )=A(index_F_global_uy_plus ,index_F_global_uy_plus )-TTrot(4,4)*(FSIe(index_F_elem_uy_plus ,index_F_elem_uy_plus ));
            
            
        otherwise
            plante_dans_ZOD_application
            
    end
end