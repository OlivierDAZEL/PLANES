% internal_fluid.m
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
e_edge=e_1;

a=nodes(edges.flux(ie,1),1:2)';
b=nodes(edges.flux(ie,2),1:2)';

%%%%% vector normal pointing out from e_1

[nx,ny]=normal_edge_out_element(a,b,c_1);

M_1=M_e;
M_2=M_e;
C_2=[nx ny 0;0 0 1];
C_1=[nx ny 0;0 0 1];

P_1_in=[nx;ny;air.Z];
Q_1_in=[nx/2 ny/2 1/(2*air.Z)];
P_1_out=[nx;ny;-air.Z];
Q_1_out=[nx/2 ny/2 -1/(2*air.Z)];

P_2_in=P_1_out;
Q_2_in=Q_1_out;
P_2_out=P_1_in;
Q_2_out=Q_1_in;

R_tilde=inv([C_1*P_1_out -C_2*P_2_out])*([-C_1*P_1_in C_2*P_2_in]);

A_x=[0 0 1/air.rho;0 0 0; air.K 0 0];

P_11_tilde=(P_1_in+P_1_out*R_tilde(1,1))*Q_1_in;
P_12_tilde=P_1_out*R_tilde(1,2)*Q_2_in;
P_21_tilde=P_2_out*R_tilde(2,1)*Q_1_in;
P_22_tilde=(P_2_in+P_2_out*R_tilde(2,2))*Q_2_in;


P_11_tilde=(P_1_in)*Q_1_in;
P_12_tilde=P_1_out*Q_2_in;
P_21_tilde=P_2_out*Q_1_in;
P_22_tilde=(P_2_in)*Q_2_in;



F_11_tilde=-(P_11_tilde(1,:)*nx+P_11_tilde(2,:)*ny)/(1j*omega);
F_12_tilde=-(P_12_tilde(1,:)*nx+P_12_tilde(2,:)*ny)/(1j*omega);
F_21_tilde=+(P_21_tilde(1,:)*nx+P_21_tilde(2,:)*ny)/(1j*omega); %+ because n2
F_22_tilde=+(P_22_tilde(1,:)*nx+P_22_tilde(2,:)*ny)/(1j*omega); %+ because n2



switch elem.model(e_1)
    case 1 % TR1
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
        vx_e1d1=-derive_polynom_2D_x(p_e1d1)/(1j*omega*rho_e);
        vx_e1d2=-derive_polynom_2D_x(p_e1d2)/(1j*omega*rho_e);
        vx_e1d3=-derive_polynom_2D_x(p_e1d3)/(1j*omega*rho_e);
        vx_e1d4=-derive_polynom_2D_x(p_e1d4)/(1j*omega*rho_e);
        vx_e1d5=-derive_polynom_2D_x(p_e1d5)/(1j*omega*rho_e);
        vx_e1d6=-derive_polynom_2D_x(p_e1d6)/(1j*omega*rho_e);
        vy_e1d1=-derive_polynom_2D_y(p_e1d1)/(1j*omega*rho_e);
        vy_e1d2=-derive_polynom_2D_y(p_e1d2)/(1j*omega*rho_e);
        vy_e1d3=-derive_polynom_2D_y(p_e1d3)/(1j*omega*rho_e);
        vy_e1d4=-derive_polynom_2D_y(p_e1d4)/(1j*omega*rho_e);
        vy_e1d5=-derive_polynom_2D_y(p_e1d5)/(1j*omega*rho_e);
        vy_e1d6=-derive_polynom_2D_y(p_e1d6)/(1j*omega*rho_e);
    case 3 % TR3
        index_p_1=dof_A(p_TR(elem.nodes(e_1,1:3)));
        nb_dof_1=3;
        vcor=nodes(nonzeros(elem.nodes(e_1,1:3)),1:2);
        p_e1d1=Lagrange_TR3(vcor,1);
        p_e1d2=Lagrange_TR3(vcor,2);
        p_e1d3=Lagrange_TR3(vcor,3);
        vx_e1d1=-derive_polynom_2D_x(p_e1d1)/(1j*omega*rho_e);
        vx_e1d2=-derive_polynom_2D_x(p_e1d2)/(1j*omega*rho_e);
        vx_e1d3=-derive_polynom_2D_x(p_e1d3)/(1j*omega*rho_e);
        vy_e1d1=-derive_polynom_2D_y(p_e1d1)/(1j*omega*rho_e);
        vy_e1d2=-derive_polynom_2D_y(p_e1d2)/(1j*omega*rho_e);
        vy_e1d3=-derive_polynom_2D_y(p_e1d3)/(1j*omega*rho_e);
    case 2
        index_p_1=dof_A(p_H12(elem.nodes(e_1,1:4)));
        nb_dof_1=12;
        vcor=nodes(nonzeros(elem.nodes(e_1,1:4)),1:2);
        base_e1=mean(vcor);       
        lx=norm(nodes(elem.nodes(e_1,1),:)-nodes(elem.nodes(e_1,2),:));
        ly=norm(nodes(elem.nodes(e_1,1),:)-nodes(elem.nodes(e_1,4),:));
        [p_e1d1,p_e1d2,p_e1d3,p_e1d4,p_e1d5,p_e1d6,p_e1d7,p_e1d8,p_e1d9,p_e1d10,p_e1d11,p_e1d12]=H12_shape_functions_shifted(lx,ly,0,0);
        
        vx_e1d1 =-derive_polynom_2D_x(p_e1d1 )/(1j*omega*rho_e);
        vx_e1d2 =-derive_polynom_2D_x(p_e1d2 )/(1j*omega*rho_e);
        vx_e1d3 =-derive_polynom_2D_x(p_e1d3 )/(1j*omega*rho_e);
        vx_e1d4 =-derive_polynom_2D_x(p_e1d4 )/(1j*omega*rho_e);
        vx_e1d5 =-derive_polynom_2D_x(p_e1d5 )/(1j*omega*rho_e);
        vx_e1d6 =-derive_polynom_2D_x(p_e1d6 )/(1j*omega*rho_e);
        vx_e1d7 =-derive_polynom_2D_x(p_e1d7 )/(1j*omega*rho_e);
        vx_e1d8 =-derive_polynom_2D_x(p_e1d8 )/(1j*omega*rho_e);
        vx_e1d9 =-derive_polynom_2D_x(p_e1d9 )/(1j*omega*rho_e);
        vx_e1d10=-derive_polynom_2D_x(p_e1d10)/(1j*omega*rho_e);
        vx_e1d11=-derive_polynom_2D_x(p_e1d11)/(1j*omega*rho_e);
        vx_e1d12=-derive_polynom_2D_x(p_e1d12)/(1j*omega*rho_e);
        
        vy_e1d1 =-derive_polynom_2D_y(p_e1d1 )/(1j*omega*rho_e);
        vy_e1d2 =-derive_polynom_2D_y(p_e1d2 )/(1j*omega*rho_e);
        vy_e1d3 =-derive_polynom_2D_y(p_e1d3 )/(1j*omega*rho_e);
        vy_e1d4 =-derive_polynom_2D_y(p_e1d4 )/(1j*omega*rho_e);
        vy_e1d5 =-derive_polynom_2D_y(p_e1d5 )/(1j*omega*rho_e);
        vy_e1d6 =-derive_polynom_2D_y(p_e1d6 )/(1j*omega*rho_e);
        vy_e1d7 =-derive_polynom_2D_y(p_e1d7 )/(1j*omega*rho_e);
        vy_e1d8 =-derive_polynom_2D_y(p_e1d8 )/(1j*omega*rho_e);
        vy_e1d9 =-derive_polynom_2D_y(p_e1d9 )/(1j*omega*rho_e);
        vy_e1d10=-derive_polynom_2D_y(p_e1d10)/(1j*omega*rho_e);
        vy_e1d11=-derive_polynom_2D_y(p_e1d11)/(1j*omega*rho_e);
        vy_e1d12=-derive_polynom_2D_y(p_e1d12)/(1j*omega*rho_e);
        
        
end
switch elem.model(e_2)
    case 1
        index_p_2=dof_A(p_TR(elem.nodes(e_2,1:6)));
        nb_dof_2=6;
        vcor=nodes(nonzeros(elem.nodes(e_2,1:6)),1:2);
        p_e2d1=Lagrange_TR6(vcor,1);
        p_e2d2=Lagrange_TR6(vcor,2);
        p_e2d3=Lagrange_TR6(vcor,3);
        p_e2d4=Lagrange_TR6(vcor,4);
        p_e2d5=Lagrange_TR6(vcor,5);
        p_e2d6=Lagrange_TR6(vcor,6);
        vx_e2d1=-derive_polynom_2D_x(p_e2d1)/(1j*omega*rho_e);
        vx_e2d2=-derive_polynom_2D_x(p_e2d2)/(1j*omega*rho_e);
        vx_e2d3=-derive_polynom_2D_x(p_e2d3)/(1j*omega*rho_e);
        vx_e2d4=-derive_polynom_2D_x(p_e2d4)/(1j*omega*rho_e);
        vx_e2d5=-derive_polynom_2D_x(p_e2d5)/(1j*omega*rho_e);
        vx_e2d6=-derive_polynom_2D_x(p_e2d6)/(1j*omega*rho_e);
        vy_e2d1=-derive_polynom_2D_y(p_e2d1)/(1j*omega*rho_e);
        vy_e2d2=-derive_polynom_2D_y(p_e2d2)/(1j*omega*rho_e);
        vy_e2d3=-derive_polynom_2D_y(p_e2d3)/(1j*omega*rho_e);
        vy_e2d4=-derive_polynom_2D_y(p_e2d4)/(1j*omega*rho_e);
        vy_e2d5=-derive_polynom_2D_y(p_e2d5)/(1j*omega*rho_e);
        vy_e2d6=-derive_polynom_2D_y(p_e2d6)/(1j*omega*rho_e);
    case 3
        index_p_2=dof_A(p_TR(elem.nodes(e_2,1:3)));
        nb_dof_2=3;
        vcor=nodes(nonzeros(elem.nodes(e_2,1:3)),1:2);
        p_e2d1=Lagrange_TR3(vcor,1);
        p_e2d2=Lagrange_TR3(vcor,2);
        p_e2d3=Lagrange_TR3(vcor,3);
        vx_e2d1=-derive_polynom_2D_x(p_e2d1)/(1j*omega*rho_e);
        vx_e2d2=-derive_polynom_2D_x(p_e2d2)/(1j*omega*rho_e);
        vx_e2d3=-derive_polynom_2D_x(p_e2d3)/(1j*omega*rho_e);
        vy_e2d1=-derive_polynom_2D_y(p_e2d1)/(1j*omega*rho_e);
        vy_e2d2=-derive_polynom_2D_y(p_e2d2)/(1j*omega*rho_e);
        vy_e2d3=-derive_polynom_2D_y(p_e2d3)/(1j*omega*rho_e);
    case 2
        index_p_2=dof_A(p_H12(elem.nodes(e_2,1:4)));
        nb_dof_2=12;
        base_e2=nodes(elem.nodes(e_2,1),:);
        lx=norm(nodes(elem.nodes(e_2,1),:)-nodes(elem.nodes(e_2,2),:));
        ly=norm(nodes(elem.nodes(e_2,1),:)-nodes(elem.nodes(e_2,4),:));
        [p_e2d1,p_e2d2,p_e2d3,p_e2d4,p_e2d5,p_e2d6,p_e2d7,p_e2d8,p_e2d9,p_e2d10,p_e2d11,p_e2d12]=H12_shape_functions(lx,ly);


        vx_e2d1 =-derive_polynom_2D_x(p_e2d1 )/(1j*omega*rho_e);
        vx_e2d2 =-derive_polynom_2D_x(p_e2d2 )/(1j*omega*rho_e);
        vx_e2d3 =-derive_polynom_2D_x(p_e2d3 )/(1j*omega*rho_e);
        vx_e2d4 =-derive_polynom_2D_x(p_e2d4 )/(1j*omega*rho_e);
        vx_e2d5 =-derive_polynom_2D_x(p_e2d5 )/(1j*omega*rho_e);
        vx_e2d6 =-derive_polynom_2D_x(p_e2d6 )/(1j*omega*rho_e);
        vx_e2d7 =-derive_polynom_2D_x(p_e2d7 )/(1j*omega*rho_e);
        vx_e2d8 =-derive_polynom_2D_x(p_e2d8 )/(1j*omega*rho_e);
        vx_e2d9 =-derive_polynom_2D_x(p_e2d9 )/(1j*omega*rho_e);
        vx_e2d10=-derive_polynom_2D_x(p_e2d10)/(1j*omega*rho_e);
        vx_e2d11=-derive_polynom_2D_x(p_e2d11)/(1j*omega*rho_e);
        vx_e2d12=-derive_polynom_2D_x(p_e2d12)/(1j*omega*rho_e);
        
        vy_e2d1 =-derive_polynom_2D_y(p_e2d1 )/(1j*omega*rho_e);
        vy_e2d2 =-derive_polynom_2D_y(p_e2d2 )/(1j*omega*rho_e);
        vy_e2d3 =-derive_polynom_2D_y(p_e2d3 )/(1j*omega*rho_e);
        vy_e2d4 =-derive_polynom_2D_y(p_e2d4 )/(1j*omega*rho_e);
        vy_e2d5 =-derive_polynom_2D_y(p_e2d5 )/(1j*omega*rho_e);
        vy_e2d6 =-derive_polynom_2D_y(p_e2d6 )/(1j*omega*rho_e);
        vy_e2d7 =-derive_polynom_2D_y(p_e2d7 )/(1j*omega*rho_e);
        vy_e2d8 =-derive_polynom_2D_y(p_e2d8 )/(1j*omega*rho_e);
        vy_e2d9 =-derive_polynom_2D_y(p_e2d9 )/(1j*omega*rho_e);
        vy_e2d10=-derive_polynom_2D_y(p_e2d10)/(1j*omega*rho_e);
        vy_e2d11=-derive_polynom_2D_y(p_e2d11)/(1j*omega*rho_e);
        vy_e2d12=-derive_polynom_2D_y(p_e2d12)/(1j*omega*rho_e);
        
        
end

for i_test=1:nb_dof_1
    eval(['Interp_test=p_e1d',num2str(i_test),';']);
    base_test=base_e1;
    for i_champs=1:nb_dof_1
        base_champs=base_e1;
        eval(['Interp_champs=F_11_tilde(1)*vx_e1d',num2str(i_champs),'+F_11_tilde(2)*vy_e1d',num2str(i_champs),'+F_11_tilde(3)*p_e1d',num2str(i_champs),';'])
        A(index_p_1(i_test),index_p_1(i_champs))=A(index_p_1(i_test),index_p_1(i_champs))+integrate_polynom_2D_edge(Interp_test,base_test,Interp_champs,base_champs,a,b,Gauss_points);
    end
    for i_champs=1:nb_dof_2
        base_champs=base_e2;
        eval(['Interp_champs=F_12_tilde(1)*vx_e2d',num2str(i_champs),'+F_12_tilde(2)*vy_e2d',num2str(i_champs),'+F_12_tilde(3)*p_e2d',num2str(i_champs),';'])
        A(index_p_1(i_test),index_p_2(i_champs))=A(index_p_1(i_test),index_p_2(i_champs))+integrate_polynom_2D_edge(Interp_test,base_test,Interp_champs,base_champs,a,b,Gauss_points);
    end
end

for i_test=1:nb_dof_2
    eval(['Interp_test=p_e2d',num2str(i_test),';']);
    base_test=base_e2;
    for i_champs=1:nb_dof_1
        base_champs=base_e1;
        eval(['Interp_champs=F_21_tilde(1)*vx_e1d',num2str(i_champs),'+F_21_tilde(2)*vy_e1d',num2str(i_champs),'+F_21_tilde(3)*p_e1d',num2str(i_champs),';'])
        A(index_p_2(i_test),index_p_1(i_champs))=A(index_p_2(i_test),index_p_1(i_champs))+integrate_polynom_2D_edge(Interp_test,base_test,Interp_champs,base_champs,a,b,Gauss_points);
    end
    for i_champs=1:nb_dof_2
        base_champs=base_e2;
        eval(['Interp_champs=F_22_tilde(1)*vx_e2d',num2str(i_champs),'+F_22_tilde(2)*vy_e2d',num2str(i_champs),'+F_22_tilde(3)*p_e2d',num2str(i_champs),';'])
        A(index_p_2(i_test),index_p_2(i_champs))=A(index_p_2(i_test),index_p_2(i_champs))+integrate_polynom_2D_edge(Interp_test,base_test,Interp_champs,base_champs,a,b,Gauss_points);
    end
end

