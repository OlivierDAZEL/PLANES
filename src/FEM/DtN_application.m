% DtN_application.m
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

nb_dof_temp=nb.dof_FEM+nb.dof_DGM;

for ie=1:nb.DtN
    typ=floor(edges.DtN(ie,4));
    
    x1=nodes(edges.DtN(ie,1),1);
    y1=nodes(edges.DtN(ie,1),2);
    x2=nodes(edges.DtN(ie,2),1);
    y2=nodes(edges.DtN(ie,2),2);
    length_edge=sqrt((x2-x1)^2+(y2-y1)^2);
    if (x1<x2)
        a=x1;
        node(1)=edges.DtN(ie,1);
        node(2)=edges.DtN(ie,2);
        node(3)=edges.DtN(ie,6);
    else
        a=x2;
        node(2)=edges.DtN(ie,1);
        node(1)=edges.DtN(ie,2);
        node(3)=edges.DtN(ie,6);
    end
    
    switch typ
        case {10}
            DtN_acou_R=1;
            F3=TR6_PW(length_edge,k_x,a);
            index_force=dof_A(p_TR(node));
            index_F_elem=find(index_force);
            index_F_global=index_force(index_F_elem);
            Omega_F.uy=1j*k_z/(air.rho*omega^2);
            F(index_F_global)=F(index_F_global)+Omega_F.uy*F3(index_F_elem);
            for i_R=1:nb.R
                F3=TR6_PW(length_edge,vec_k_x(i_R),a);
                index_force=dof_A(p_TR(node));
                index_F_elem=find(index_force);
                index_F_global=index_force(index_F_elem);
                A(index_F_global,nb_dof_temp+i_R)=A(index_F_global,nb_dof_temp+i_R)+F3(index_F_elem)*(1j*vec_k_z(i_R)/(air.rho*omega^2));
                A(nb_dof_temp+i_R,index_F_global)=A(nb_dof_temp+i_R,index_F_global)+F3(index_F_elem)';
            end
        case {11}
            DtN_elas_R=1;
            F3=TR6_PW(length_edge,k_x,a);
            index_force=dof_A(uy_TR(node));
            index_F_elem=find(index_force);
            index_F_global=index_force(index_F_elem);
            Omega_F.p=1;
            F(index_F_global)=F(index_F_global)+Omega_F.p*F3(index_F_elem);
            for i_R=1:nb.R
                F3=TR6_PW(length_edge,vec_k_x(i_R),a);
                index_force=dof_A(uy_TR(node));
                index_F_elem=find(index_force);
                index_F_global=index_force(index_F_elem);
                A(index_F_global,nb_dof_temp+i_R)=A(index_F_global,nb_dof_temp+i_R)-F3(index_F_elem);
                A(nb_dof_temp+i_R,index_F_global)=A(nb_dof_temp+i_R,index_F_global)+F3(index_F_elem)';
            end
        case {13}
            DtN_plate_R=1;
            Omega_moins=[0;-1j*k_z/(air.rho*omega^2);-1;0];
            Omega_plus=transfer_unknowns(k_x,omega,Omega_moins,1,data_model.incident);
            F3=TR6_PW(length_edge,k_x,a);
            index_force=dof_A(ux_TR(node));
            index_F_elem=find(index_force);
            index_F_global=index_force(index_F_elem);
            F(index_F_global)=F(index_F_global)-Omega_plus(1,1)*F3(index_F_elem);
            
            index_force=dof_A(uy_TR(node));
            index_F_elem=find(index_force);
            index_F_global=index_force(index_F_elem);
            F(index_F_global)=F(index_F_global)-Omega_plus(3,1)*F3(index_F_elem);
            
            for i_R=1:nb.R
                
                temp=nb_dof_temp+1+size_info_vector_R*(i_R-1);
                
                F3=TR6_PW(length_edge,vec_k_x(i_R),a);
                
                Omega_moins=[0 0;1j*vec_k_z(i_R)/(air.rho*omega^2) 0;-1 0;0 1];
                [Omega_plus]=transfer_unknowns(vec_k_x(i_R),omega,Omega_moins,1,data_model.incident);
                

                index_force=dof_A(ux_TR(node));
                index_F_elem=find(index_force);
                index_F_global=index_force(index_F_elem);
                
                A(index_F_global,temp  )=A(index_F_global,temp  )+Omega_plus(1,1)*F3(index_F_elem);
                A(index_F_global,temp+1)=A(index_F_global,temp+1)+Omega_plus(1,2)*F3(index_F_elem);
                A(temp  ,index_F_global)=A(temp,index_F_global)+F3(index_F_elem)';
                
                index_force=dof_A(uy_TR(node));
                index_F_elem=find(index_force);
                index_F_global=index_force(index_F_elem);
                A(index_F_global,temp  )=A(index_F_global,temp  )+Omega_plus(3,1)*F3(index_F_elem);
                A(index_F_global,temp+1)=A(index_F_global,temp+1)+Omega_plus(3,2)*F3(index_F_elem);
                A(temp+1,index_F_global)=A(temp+1,index_F_global)+F3(index_F_elem)';
            end
            
        case {12}
            DtN_2001_R=1;
            F3=TR6_PW(length_edge,k_x,a);
            
            % Terme p_a delta u_y champs incident
            index_force=dof_A(uy_TR(node));
            index_F_elem=find(index_force);
            index_F_global=index_force(index_F_elem);
            F(index_F_global)=F(index_F_global)+F3(index_F_elem);
            %signe ok
            
            % Terme u_a delta p champs incident
            index_force=dof_A(p_TR(node));
            index_F_elem=find(index_force);
            index_F_global=index_force(index_F_elem);
            F(index_F_global)=F(index_F_global)+F3(index_F_elem)*(1i*k_z)/(air.rho*omega^2);
            %signe ok
            
            %! Reflected fields
            
            for i_R=1:nb.R
                % Terme p_a delta u_y champs reflechi
                
                F3=TR6_PW(length_edge,vec_k_x(i_R),a);
                index_force=dof_A(uy_TR(node));
                index_F_elem=find(index_force);
                index_F_global=index_force(index_F_elem);
                A(index_F_global,nb_dof_temp+i_R)=A(index_F_global,nb_dof_temp+i_R)-F3(index_F_elem);
                
                
                % Terme u_a delta p champs reflechi
                
                index_force=dof_A(p_TR(node));
                index_F_elem=find(index_force);
                index_F_global=index_force(index_F_elem);
                A(index_F_global,nb_dof_temp+i_R)=A(index_F_global,nb_dof_temp+i_R)+F3(index_F_elem)*(1i*vec_k_z(i_R))/(air.rho*omega^2);
                
                %%%%%%% Equation suppl?mentare sur la pression
                A(nb_dof_temp+i_R,index_F_global)=A(nb_dof_temp+i_R,index_F_global)+F3(index_F_elem)';
            end
                        
            a1(1)=nodes(node(1),1);
            a1(2)=nodes(node(1),2);
            a2(1)=nodes(node(2),1);
            a2(2)=nodes(node(2),2);
            
            FSIe=TR6_FSI(a1,a2);
            
            
            index_force_p=dof_A(p_TR(node));
            index_F_elem_p=find(index_force_p);
            index_F_global_p=index_force_p(index_F_elem_p);
            
            index_force_u=dof_A(uy_TR(node));
            index_F_elem_u=find(index_force_u);
            index_F_global_u=index_force_u(index_F_elem_u);
            
            A(index_F_global_p,index_F_global_u)=A(index_F_global_p,index_F_global_u)-(FSIe(index_F_elem_p,index_F_elem_u));
        case {20}
            DtN_acou_T=1;
            for i_T=1:nb.T
                F3=TR6_PW(length_edge,vec_k_x_t(i_T),a);
                index_force=dof_A(p_TR(node));
                index_F_elem=find(index_force);
                index_F_global=index_force(index_F_elem);
                Omega.uy=1j*vec_k_z_t(i_T)/(air.rho*omega^2);
                A(index_F_global,nb_dof_temp+size_info_vector_R*nb.R+i_T)=A(index_F_global,nb_dof_temp+size_info_vector_R*nb.R+i_T)+Omega.uy*F3(index_F_elem);
                A(nb_dof_temp+size_info_vector_R*nb.R+i_T,index_F_global)=A(nb_dof_temp+size_info_vector_R*nb.R+i_T,index_F_global)+         F3(index_F_elem)';
            end
            
        case {21}
            DtN_elas_T=1;
            F3=TR6_PW(length_edge,k_x,a);
            for i_T=1:nb.T
                F3=TR6_PW(length_edge,vec_k_x(i_T),a);
                index_force=dof_A(uy_TR(node));
                index_F_elem=find(index_force);
                index_F_global=index_force(index_F_elem);
                Omega.p=1;
                A(index_F_global,nb_dof_temp+size_info_vector_R*nb.R+i_T)=A(index_F_global,nb_dof_temp+size_info_vector_R*nb.R+i_T)+Omega.p*F3(index_F_elem);
                A(nb_dof_temp+size_info_vector_R*nb.R+i_T,index_F_global)=A(nb_dof_temp+size_info_vector_R*nb.R+i_T,index_F_global)+        F3(index_F_elem)';
            end
        case {23}
            
            DtN_plate_T=1;
            
            for i_T=1:nb.T
                temp=nb_dof_temp+size_info_vector_R*nb.R+1+size_info_vector_T*(i_T-1);
                
                F3=TR6_PW(length_edge,vec_k_x_t(i_T),a);
                
                Omega_moins=[0 0;-1j*vec_k_z(i_T)/(air.rho*omega^2) 0;-1 0;0 1];
                [Omega_plus]=transfer_unknowns(vec_k_x_t(i_T),omega,Omega_moins,-1,data_model.transmitted);
                
                index_force=dof_A(ux_TR(node));
                index_F_elem=find(index_force);
                index_F_global=index_force(index_F_elem);
                A(index_F_global,temp  )=A(index_F_global,temp  )-Omega_plus(1,1)*F3(index_F_elem);
                A(index_F_global,temp+1)=A(index_F_global,temp+1)-Omega_plus(1,2)*F3(index_F_elem);
                A(temp,index_F_global)=A(temp,index_F_global)+F3(index_F_elem)';
                
                index_force=dof_A(uy_TR(node));
                index_F_elem=find(index_force);
                index_F_global=index_force(index_F_elem);
                
                A(index_F_global,temp  )=A(index_F_global,temp  )-Omega_plus(3,1)*F3(index_F_elem);
                A(index_F_global,temp+1)=A(index_F_global,temp+1)-Omega_plus(3,2)*F3(index_F_elem);
                
                A(temp+1,index_F_global)=A(temp+1,index_F_global)+F3(index_F_elem)';
                
            end
            
            
        otherwise
            disp('Unknown load')
            stop
    end
end

if exist('DtN_acou_R')
    for i_R=1:nb.R
        A(i_R+nb_dof_temp,i_R+nb_dof_temp)=-period;
    end
    F(nb_dof_temp+1)=F(nb_dof_temp+1)+period;
end

if exist('DtN_elas_R')
    for i_R=1:nb.R
        A(i_R+nb_dof_temp,i_R+nb_dof_temp)=-period*(1i*vec_k_z(i_R))/(air.rho*omega^2);
    end
    F(nb_dof_temp+1)=F(nb_dof_temp+1)-period*(1i*k_z)/(air.rho*omega^2);
end

if exist('DtN_plate_R')
    
    Omega_moins=[0;-1j*k_z/(air.rho*omega^2);-1 ;0 ];
    Omega_plus=transfer_unknowns(k_x,omega,Omega_moins,1,data_model.incident);
    
    F(nb_dof_temp+1)=F(nb_dof_temp+1)+period*Omega_plus(4);
    F(nb_dof_temp+2)=F(nb_dof_temp+2)+period*Omega_plus(2);
    
    for i_R=1:nb.R
        temp=nb_dof_temp+1+size_info_vector_R*(i_R-1);    
        Omega_moins=[0 0;1j*vec_k_z(i_R)/(air.rho*omega^2) 0;-1 0;0 1];
        [Omega_plus]=transfer_unknowns(vec_k_x(i_R),omega,Omega_moins,1,data_model.incident);
        A(temp  ,temp  )=-period*Omega_plus(4,1);
        A(temp  ,temp+1)=-period*Omega_plus(4,2);
        A(temp+1,temp  )=-period*Omega_plus(2,1);
        A(temp+1,temp+1)=-period*Omega_plus(2,2);
    end
end

if exist('DtN_plate_T')
    
    for i_T=1:nb.T
        temp=nb_dof_temp+size_info_vector_R*nb.R+1+size_info_vector_T*(i_T-1);
        Omega_moins=[0 0;-1j*vec_k_z_t(i_T)/(air.rho*omega^2) 0;-1 0;0 1];
        [Omega_plus]=transfer_unknowns(vec_k_x_t(i_T),omega,Omega_moins,-1,data_model.transmitted);
        A(temp  ,temp  )=-period*Omega_plus(4,1);
        A(temp  ,temp+1)=-period*Omega_plus(4,2);
        A(temp+1,temp  )=-period*Omega_plus(2,1);
        A(temp+1,temp+1)=-period*Omega_plus(2,2);
        
    end
end

if exist('DtN_acou_T')
    for i_T=1:nb.T
        A(size_info_vector_R*nb.R+i_T+nb_dof_temp,size_info_vector_R*nb.R+i_T+nb_dof_temp)=-period;
    end
end

if exist('DtN_elas_T')
    for i_T=1:nb.T
        A(size_info_vector_R*nb.R+i_T+nb_dof_temp,size_info_vector_R*nb.R+i_T+nb_dof_temp)=period*(1i*vec_k_z_t(i_T))/(air.rho*omega^2);
    end
end

if exist('DtN_2001_R')
    F(nb_dof_temp+1)=period;
    for i_R=1:nb.R
        A(nb_dof_temp+i_R,nb_dof_temp+i_R)=A(nb_dof_temp+i_R,nb_dof_temp+i_R)-period;
    end
end







