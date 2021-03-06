% internal_fluid_FEM_DGM.m
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
%                 https://github.com/OlivierDAZEL/PLANES
% or find more details on Olivier's webpage
%                  http://perso.univ-lemans.fr/~odazel/
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
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


c_e=air.c;
k_e=omega/c_e;
Z_e=air.Z;
rho_e=air.rho;
M_e=diag([air.rho,air.rho,1/air.K]);


e_1=edges.flux(ie,3);
e_2=edges.flux(ie,4);
c_1=mean(nodes(nonzeros(elem.nodes(e_1,:)),1:2))';
c_2=mean(nodes(nonzeros(elem.nodes(e_2,:)),1:2))';

a=nodes(edges.flux(ie,1),1:2)';
b=nodes(edges.flux(ie,2),1:2)';

%%%%% vector normal pointing out from e_1
[nx,ny]=normal_edge_out_element(a,b,c_1);

M_2=M_e;
C_2=[nx ny 0;0 0 1];

P_1_in=[nx;ny;air.Z];
Q_1_in=[nx/2 ny/2 1/(2*air.Z)];
P_1_out=[nx;ny;-air.Z];
Q_1_out=[nx/2 ny/2 -1/(2*air.Z)];

P_2_in=P_1_out;
Q_2_in=Q_1_out;
P_2_out=P_1_in;
Q_2_out=Q_1_in;

R_tilde=inv([[1; 0], -C_2*P_2_out])*([[0; -1], C_2*P_2_in]);

A_x=[0 0 1/air.rho;0 0 0; air.K 0 0];
A_y=[0 0 0;0 0 1/air.rho; 0 air.K 0];
F_2=-(A_x*nx+A_y*ny);

F_11_tilde=-R_tilde(1,1);%/(1j*omega);
F_12_tilde=-R_tilde(1,2)*Q_2_in;%/(1j*omega);

F_21_tilde=M_2*F_2*P_2_out*R_tilde(2,1);
F_22_tilde=M_2*F_2*(P_2_in + P_2_out*R_tilde(2,2))*Q_2_in;


switch elem.model(e_1)
    case 1 % TR6
        index_p_1=dof_A(p_TR(elem.nodes(e_1,1:6)));
        nb_dof_1=6;
        vcor=nodes(nonzeros(elem.nodes(e_1,1:6)),1:2);
        base_e1=mean(vcor);
        p_e1d1=Lagrange_TR6(vcor,1,base_e1);
        p_e1d2=Lagrange_TR6(vcor,2,base_e1);
        p_e1d3=Lagrange_TR6(vcor,3,base_e1);
        p_e1d4=Lagrange_TR6(vcor,4,base_e1);
        p_e1d5=Lagrange_TR6(vcor,5,base_e1);
        p_e1d6=Lagrange_TR6(vcor,6,base_e1);
    case 3 % TR3
        index_p_1=dof_A(p_TR(elem.nodes(e_1,1:3)));
        nb_dof_1=3;
        vcor=nodes(nonzeros(elem.nodes(e_1,1:3)),1:2);
        p_e1d1=Lagrange_TR3(vcor,1);
        p_e1d2=Lagrange_TR3(vcor,2);
        p_e1d3=Lagrange_TR3(vcor,3);
    case 2
        index_p_1=dof_A(p_H12(elem.nodes(e_1,1:4)));
        nb_dof_1=12;
        vcor=nodes(nonzeros(elem.nodes(e_1,1:4)),1:2);
        base_e1=nodes(elem.nodes(e_1,1),:);
        lx=norm(nodes(elem.nodes(e_1,1),:)-nodes(elem.nodes(e_1,2),:));
        ly=norm(nodes(elem.nodes(e_1,1),:)-nodes(elem.nodes(e_1,4),:));
        [p_e1d1,p_e1d2,p_e1d3,p_e1d4,p_e1d5,p_e1d6,p_e1d7,p_e1d8,p_e1d9,p_e1d10,p_e1d11,p_e1d12]=H12_shape_functions_shifted(lx,ly,0,0);
end

index_2  =((1:data_model.theta_DGM.nb)-1)+dof_start_element(e_2);

for i_test=1:nb_dof_1
    eval(['N_1=p_e1d',num2str(i_test),';']);
    base_test=base_e1;
    for i_champs=1:nb_dof_1
        base_champs=base_e1;
        eval(['I_tilde = F_11_tilde*p_e1d', num2str(i_champs), ';']);
        A(index_p_1(i_test),index_p_1(i_champs))=A(index_p_1(i_test),index_p_1(i_champs))+integrate_polynom_2D_edge(N_1,base_test,I_tilde,base_champs,a,b,Gauss_points);
    end
    for i_champs=1:data_model.theta_DGM.nb
        jk=-1j*k_e*[cos(vec_theta(i_champs));sin(vec_theta(i_champs))];
        Phi_champs=Phi_fluid(cos(vec_theta(i_champs)),sin(vec_theta(i_champs)),Z_e);
        I_tilde=F_12_tilde*Phi_champs;
        A(index_p_1(i_test),index_2(i_champs))=A(index_p_1(i_test),index_2(i_champs))+I_tilde*integrate_polynom_exp_2D_edge(N_1,base_test,jk,c_2,a,b,Gauss_points);
    end
end


for i_test=1:data_model.theta_DGM.nb
    Phi_test=Phi_fluid(cos(vec_theta(i_test)),sin(vec_theta(i_test)),Z_e);
    for i_champs=1:nb_dof_1
        base_champs=base_e1;
        delta_exp=exp(1j*k_e*[cos(vec_theta(i_test));sin(vec_theta(i_test))]'*(a-c_2));
        jk=1j*k_e*[cos(vec_theta(i_test));sin(vec_theta(i_test))];

        eval(['N_1 = p_e1d', num2str(i_champs), ';']);
        I_tilde=Phi_test.'*F_21_tilde;
        A(index_2(i_test),index_p_1(i_champs))=A(index_2(i_test),index_p_1(i_champs))+I_tilde*integrate_polynom_exp_2D_edge(N_1,base_champs,jk,c_2,a,b,Gauss_points);
    end
end

nx=cos(vec_theta);
ny=sin(vec_theta);
Phi=Phi_fluid_vector(nx,ny,Z_e,Shift_fluid);
II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a,b,[c_2 c_2]);
MM=kron(II,F_22_tilde);
A(index_2,index_2)=A(index_2,index_2)+Phi.'*MM*Phi;
