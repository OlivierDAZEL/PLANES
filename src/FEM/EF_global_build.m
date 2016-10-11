% EF_TR6_global_build.m
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

create_temporary_FEM_matrices

project.logger(1, 'FEM', 'Beginning construction of FE related matrices');

tic
for ie=1:nb.elements
      
    if elem.model(ie)==1 % TR6
        nodes_elements = nodes(elem.nodes(ie,1:6),1:2)';
        if floor(elem.label(ie)/1000)==0
            [vh,vq]=TR6_fluid(nodes_elements);
            index_p=p_TR(elem.nodes(ie,1:6));
            H_acou(index_p,index_p)=H_acou(index_p,index_p)+vh;
            Q_acou(index_p,index_p)=Q_acou(index_p,index_p)+vq;
        elseif floor(elem.label(ie)/1000)==1
            [vm,vk0,vk1]=TR6_elas(nodes_elements);
            index_e=uxy_TR(elem.nodes(ie,1:6));
            insert_temporary_matrices_elas
        elseif floor(elem.label(ie)/1000)==2
            [vh,vq]=TR6_fluid(nodes_elements);
            index_p=p_TR(elem.nodes(ie,1:6));
            insert_temporary_matrices_eqf
        elseif floor(elem.label(ie)/1000)==3
            [vh,vq]=TR6_fluid(nodes_elements);
            index_p=p_TR(elem.nodes(ie,1:6));
            insert_temporary_matrices_limp
        elseif floor(elem.label(ie)/1000)==4
            [vm,vk0,vk1,vc,vh,vq]=TR6_pem98(nodes_elements);
            index_e=uxy_TR(elem.nodes(ie,1:6));
            index_p=p_TR(elem.nodes(ie,1:6));
            insert_temporary_matrices_pem1998
        elseif floor(elem.label(ie)/1000)==5
            [vm,vk0,vk1,vc,vcp,vh,vq]=TR6_pem01(nodes_elements);
            index_e=uxy_TR(elem.nodes(ie,1:6));
            index_p=p_TR(elem.nodes(ie,1:6));
            insert_temporary_matrices_pem2001
        elseif floor(elem.label(ie)/1000)==8
%             temp=elem.label(ie);
%             temp_x=floor((temp-8000)/100);
%             temp_y=floor((temp-8000-100*temp_x)/10);
%             if temp_x==0;
%                 tau_x=1;
%             else
%                 tau_x=exp(1j*pi/4);
%             end
%             if temp_y==0;
%                 tau_y=1;
%             else
%                 tau_y=exp(1j*pi/4);
%             end
%             [vh,vq]=TR6_PML(nodes_elements,tau_x,tau_y);
%             index_p=p_TR(elem.nodes(ie,1:6));
%             H_PML(index_p,index_p)=H_PML(index_p,index_p)+vh;
%             Q_PML(index_p,index_p)=Q_PML(index_p,index_p)+vq;
        else
            disp('Subroutine parameter_element')
            disp('Unknwon fluid type of element')
            stop
        end
    elseif elem.model(ie)==3 % TR3
        nodes_elements = nodes(elem.nodes(ie,1:3),1:2)';
        if floor(elem.label(ie)/1000)==0
            [vh,vq]=TR3_fluid(nodes_elements);
            index_p=p_TR(elem.nodes(ie,1:3));
            H_acou(index_p,index_p)=H_acou(index_p,index_p)+vh;
            Q_acou(index_p,index_p)=Q_acou(index_p,index_p)+vq;
        elseif floor(elem.label(ie)/1000)==1
            [vm,vk0,vk1]=TR3_elas(nodes_elements);
            index_e=uxy_TR(elements(ie,1:6));
            insert_temporary_matrices_elas
        elseif floor(elem.label(ie)/1000)==2
            [vh,vq]=TR3_fluid(nodes_elements);
            index_p=p_TR(elements(ie,1:6));
            insert_temporary_matrices_eqf
        elseif floor(elem.label(ie)/1000)==3
            [vh,vq]=TR3_fluid(nodes_elements);
            index_p=p_TR(elements(ie,1:6));
            insert_temporary_matrices_limp
        elseif floor(elem.label(ie)/1000)==4
            [vm,vk0,vk1,vc,vh,vq]=TR3_pem98(nodes_elements);
            index_e=uxy_TR(elements(ie,1:6));
            index_p=p_TR(elements(ie,1:6));
            insert_temporary_matrices_pem1998
        elseif floor(elem.label(ie)/1000)==5
            [vm,vk0,vk1,vc,vcp,vh,vq]=TR3_pem01(nodes_elements);
            index_e=uxy_TR(elements(ie,1:6));
            index_p=p_TR(elements(ie,1:6));
            insert_temporary_matrices_pem2001
        elseif floor(elem.label(ie)/1000)==8
            temp=element_label(ie);
            temp_x=floor((temp-8000)/100);
            temp_y=floor((temp-8000-100*temp_x)/10);
            if temp_x==0;
                tau_x=1;
            else
                tau_x=exp(1j*pi/4);
            end
            if temp_y==0;
                tau_y=1;
            else
                tau_y=exp(1j*pi/4);
            end
            [vh,vq]=TR6_PML(nodes_elements,tau_x,tau_y);
            index_p=p_TR(elements(ie,1:6));
            H_PML(index_p,index_p)=H_PML(index_p,index_p)+vh;
            Q_PML(index_p,index_p)=Q_PML(index_p,index_p)+vq;
        else
            disp('Subroutine parameter_element')
            disp('Unknwon fluid type of element')
            stop
        end
        elseif elem.model(ie)==2 % H12
        nodes_elements = nodes(elem.nodes(ie,1:3),1:2)';
        if floor(elem.label(ie)/1000)==0
            lx=norm(nodes(elem.nodes(ie,1),:)-nodes(elem.nodes(ie,2),:));
            ly=norm(nodes(elem.nodes(ie,1),:)-nodes(elem.nodes(ie,4),:));

            
            %[vh,vq]=H12_fluid(lx,ly);
            index_p=p_H12(elem.nodes(ie,1:4));
            H_acou(index_p,index_p)=H_acou(index_p,index_p)+H_elem_H12(:,:,elem.H12(ie));
            Q_acou(index_p,index_p)=Q_acou(index_p,index_p)+Q_elem_H12(:,:,elem.H12(ie));
        end
    elseif elem.model(ie)==4 % TR6 axi        
        nodes_elements = nodes(elem.nodes(ie,1:6),1:2)';
        if floor(elem.label(ie)/1000)==0
            [vh,vq]=TR6_fluid_axi(nodes_elements);
            index_p=p_TR(elem.nodes(ie,1:6));
            H_acou(index_p,index_p)=H_acou(index_p,index_p)+vh;
            Q_acou(index_p,index_p)=Q_acou(index_p,index_p)+vq;
        elseif floor(elem.label(ie)/1000)==8
%             temp=elem.label(ie);
%             temp_x=floor((temp-8000)/100);
%             temp_y=floor((temp-8000-100*temp_x)/10);
%             if temp_x==0;
%                 tau_x=1;
%             else
%                 tau_x=exp(1j*pi/4);
%             end
%             if temp_y==0;
%                 tau_y=1;
%             else
%                 tau_y=exp(1j*pi/4);
%             end
%             [vh,vq]=TR6_PML_axi(nodes_elements,tau_x,tau_y);
%             index_p=p_TR(elem.nodes(ie,1:6));
%             H_PML(index_p,index_p)=H_PML(index_p,index_p)+vh;
%             Q_PML(index_p,index_p)=Q_PML(index_p,index_p)+vq;
        else
            disp('Subroutine parameter_element')
            disp('Unknwon fluid type of element')
            stop
        end
        
        
        
    end   %if on elem.model(ie)
    
    
    
    
end

size_global_matrices=12*nb.nodes;

discard_l1_temporary_FEM_matrices

H_acou=H_acou(list_dof_valid,list_dof_valid);
Q_acou=Q_acou(list_dof_valid,list_dof_valid);

% H_PML=H_PML(list_dof_valid,list_dof_valid);
% Q_PML=Q_PML(list_dof_valid,list_dof_valid);

for i_mat=1:nb.media.elas
    
    eval(['K0_elas_',num2str(i_mat),'=sparse( i_k0_elas,j_k0_elas,v_k0_elas(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    eval(['K1_elas_',num2str(i_mat),'=sparse( i_k1_elas,j_k1_elas,v_k1_elas(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    eval(['M_elas_',num2str(i_mat),'=sparse( i_m_elas,  j_m_elas,  v_m_elas(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    
    eval(['K0_elas_',num2str(i_mat),'=K0_elas_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    eval(['K1_elas_',num2str(i_mat),'=K1_elas_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    eval(['M_elas_',num2str(i_mat),'=M_elas_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    
end

for i_mat=1:nb.media.eqf
    
    eval(['Q_eqf_',num2str(i_mat),'=sparse( i_q_eqf,j_q_eqf,v_q_eqf(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    eval(['H_eqf_',num2str(i_mat),'=sparse( i_h_eqf,j_h_eqf,v_h_eqf(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    
    eval(['Q_eqf_',num2str(i_mat),'=Q_eqf_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    eval(['H_eqf_',num2str(i_mat),'=H_eqf_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    
end

for i_mat=1:nb.media.limp
    
    eval(['Q_limp_',num2str(i_mat),'=sparse( i_q_limp,j_q_limp,v_q_limp(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    eval(['H_limp_',num2str(i_mat),'=sparse( i_h_limp,j_h_limp,v_h_limp(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    
    eval(['Q_limp_',num2str(i_mat),'=Q_limp_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    eval(['H_limp_',num2str(i_mat),'=H_limp_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    
end


for i_mat=1:nb.media.pem98
    
    eval(['K0_pem98_',num2str(i_mat),'=sparse( i_k0_pem98,j_k0_pem98,v_k0_pem98(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    eval(['K1_pem98_',num2str(i_mat),'=sparse( i_k1_pem98,j_k1_pem98,v_k1_pem98(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    eval(['M_pem98_',num2str(i_mat),'=sparse( i_m_pem98,  j_m_pem98,  v_m_pem98(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    eval(['H_pem98_',num2str(i_mat),' =sparse( i_h_pem98,  j_h_pem98,  v_h_pem98(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    eval(['Q_pem98_',num2str(i_mat),' =sparse( i_q_pem98,  j_q_pem98,  v_q_pem98(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    eval(['C_pem98_',num2str(i_mat),' =sparse( i_c_pem98,  j_c_pem98,  v_c_pem98(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    
    eval(['K0_pem98_',num2str(i_mat),'=K0_pem98_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    eval(['K1_pem98_',num2str(i_mat),'=K1_pem98_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    eval(['M_pem98_',num2str(i_mat),'=M_pem98_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    eval(['H_pem98_',num2str(i_mat),'=H_pem98_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    eval(['Q_pem98_',num2str(i_mat),'=Q_pem98_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    eval(['C_pem98_',num2str(i_mat),'=C_pem98_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    
    
end


for i_mat=1:nb.media.pem01
    
    eval(['K0_pem01_',num2str(i_mat),'=sparse( i_k0_pem01,j_k0_pem01,v_k0_pem01(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    eval(['K1_pem01_',num2str(i_mat),'=sparse( i_k1_pem01,j_k1_pem01,v_k1_pem01(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    eval(['M_pem01_',num2str(i_mat),'=sparse( i_m_pem01,  j_m_pem01,  v_m_pem01(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    eval(['H_pem01_',num2str(i_mat),' =sparse( i_h_pem01,  j_h_pem01,  v_h_pem01(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    eval(['Q_pem01_',num2str(i_mat),' =sparse( i_q_pem01,  j_q_pem01,  v_q_pem01(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    eval(['C_pem01_',num2str(i_mat),' =sparse( i_c_pem01,  j_c_pem01,  v_c_pem01(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    eval(['Cp_pem01_',num2str(i_mat),'=sparse(i_cp_pem01,j_cp_pem01,   v_cp_pem01(:,',num2str(i_mat),'),size_global_matrices,size_global_matrices);']);
    
    eval(['K0_pem01_',num2str(i_mat),'=K0_pem01_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    eval(['K1_pem01_',num2str(i_mat),'=K1_pem01_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    eval(['M_pem01_',num2str(i_mat),'=M_pem01_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    eval(['H_pem01_',num2str(i_mat),'=H_pem01_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    eval(['Q_pem01_',num2str(i_mat),'=Q_pem01_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    eval(['C_pem01_',num2str(i_mat),'=C_pem01_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    eval(['Cp_pem01_',num2str(i_mat),'=Cp_pem01_',num2str(i_mat),'(list_dof_valid,list_dof_valid);']);
    
    
end


clear_temporary_FEM_matrices

_etime = toc;
project.logger(2, 'profiling', ['EF_global_build ' num2str(_etime) 's.']);
