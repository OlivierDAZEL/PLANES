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

e_1=edges.internal(ie,3);
e_2=edges.internal(ie,4);
c_1=mean(nodes(nonzeros(elem.nodes(e_1,:)),1:2))';
c_2=mean(nodes(nonzeros(elem.nodes(e_2,:)),1:2))';
e_edge=e_1;


a=nodes(edges.internal(ie,1),1:2)';
b=nodes(edges.internal(ie,2),1:2)';

%%%%% vector normal pointing out from e_1

centre_edge=(a+b)/2;
n_centre=c_1-centre_edge;
[nx,ny]=normal_edge_out_element(a,b,c_1);


% % u=[vx vy p]
% %i*omega*u+A_x du/dx+A_y dudy
%
% A_x=[0 1/air.rho 0;0 0 0; air.K 0 0];
% A_y=[0 0 0;0 0 1/air.rho; 0 air.K 0];
% F_e=(A_x*nx+A_y*ny)
% %
% C_2=[nx ny 0;0 0 1];
% C_1=[nx ny 0;0 0 1];

% W_e2=[nx -ny -nx;ny nx -ny;air.Z 0 air.Z];
%
%
% Omega_e2=[nx/2 ny/2 1/(2*air.Z);-ny nx 0;-nx/2 -ny/2 1/(2*air.Z)];
% W_e2_plus=W_e2(:,1);
% W_e2_moins=W_e2(:,3);
% Omega_e2_plus=Omega_e2(1,:);
% Omega_e2_moins=Omega_e2(3,:);
%
%
% air.c*(W_e2_plus*Omega_e2_plus-W_e2_moins*Omega_e2_moins)

W_e1_in=[nx;ny;air.Z];
Omega_e1_in=[nx/2 ny/2 1/(2*air.Z)];
W_e1_out=[nx;ny;-air.Z];
Omega_e1_out=[nx/2 ny/2 -1/(2*air.Z)];

W_e2_in=W_e1_out;
Omega_e2_in=Omega_e1_out;
W_e2_out=W_e1_in;
Omega_e2_out=Omega_e1_in;

% air.c*(W_e1_in*Omega_e1_in-W_e1_out*Omega_e1_out)
% W_e1_in*Omega_e1_in+W_e1_out*Omega_e1_out
% fsdfsdfdsfdsdffds
%inv([C_1*W_e1_out -C_2*W_e2_out])*[-C_1*W_e1_in C_2*W_e2_in]


switch elem.model(e_1)
    case 1 % TR1
        index_p_1=dof_A(p_TR(elem.nodes(e_1,1:6)));
        nb_dof_1=6;
        vcor=nodes(nonzeros(elem.nodes(e_1,1:6)),1:2)
        p_e1d1=Lagrange_TR6(vcor,1);
        p_e1d2=Lagrange_TR6(vcor,2);
        p_e1d3=Lagrange_TR6(vcor,3);
        p_e1d4=Lagrange_TR6(vcor,4);
        p_e1d5=Lagrange_TR6(vcor,5);
        p_e1d6=Lagrange_TR6(vcor,6);
        vx_e1d1=-derive_polynom_2D_x_2(p_e1d1)/(1j*omega*rho_e);
        vx_e1d2=-derive_polynom_2D_x_2(p_e1d2)/(1j*omega*rho_e);
        vx_e1d3=-derive_polynom_2D_x_2(p_e1d3)/(1j*omega*rho_e);
        vx_e1d4=-derive_polynom_2D_x_2(p_e1d4)/(1j*omega*rho_e);
        vx_e1d5=-derive_polynom_2D_x_2(p_e1d5)/(1j*omega*rho_e);
        vx_e1d6=-derive_polynom_2D_x_2(p_e1d6)/(1j*omega*rho_e);
        vy_e1d1=-derive_polynom_2D_y_2(p_e1d1)/(1j*omega*rho_e);
        vy_e1d2=-derive_polynom_2D_y_2(p_e1d2)/(1j*omega*rho_e);
        vy_e1d3=-derive_polynom_2D_y_2(p_e1d3)/(1j*omega*rho_e);
        vy_e1d4=-derive_polynom_2D_y_2(p_e1d4)/(1j*omega*rho_e);
        vy_e1d5=-derive_polynom_2D_y_2(p_e1d5)/(1j*omega*rho_e);
        vy_e1d6=-derive_polynom_2D_y_2(p_e1d6)/(1j*omega*rho_e);
    case 3 % TR3
        index_p_1=dof_A(p_TR(elem.nodes(e_1,1:3)));
        nb_dof_1=3;
        vcor=nodes(nonzeros(elem.nodes(e_1,1:3)),1:2);
        p_e1d1=Lagrange_TR3(vcor,1);
        p_e1d2=Lagrange_TR3(vcor,2);
        p_e1d3=Lagrange_TR3(vcor,3);
        vx_e1d1=-derive_polynom_2D_x_2(p_e1d1)/(1j*omega*rho_e);
        vx_e1d2=-derive_polynom_2D_x_2(p_e1d2)/(1j*omega*rho_e);
        vx_e1d3=-derive_polynom_2D_x_2(p_e1d3)/(1j*omega*rho_e);
        vy_e1d1=-derive_polynom_2D_y_2(p_e1d1)/(1j*omega*rho_e);
        vy_e1d2=-derive_polynom_2D_y_2(p_e1d2)/(1j*omega*rho_e);
        vy_e1d3=-derive_polynom_2D_y_2(p_e1d3)/(1j*omega*rho_e);
    case 2
        index_p_1=dof_A(p_H12(elem.nodes(e_1,1:4)));
        nb_dof_1=12;
        
        lx=norm(nodes(elem.nodes(e_1,1),:)-nodes(elem.nodes(e_1,2),:));
        ly=norm(nodes(elem.nodes(e_1,1),:)-nodes(elem.nodes(e_1,4),:));
        [p_e1d1,p_e1d2,p_e1d3,p_e1d4,p_e1d5,p_e1d6,p_e1d7,p_e1d8,p_e1d9,p_e1d10,p_e1d11,p_e1d12]=H12_shape_functions_shifted(lx,ly,nodes(elem.nodes(e_1,1),1),nodes(elem.nodes(e_1,1),2));
        
        vx_e1d1 =-derive_polynom_2D_x_2(p_e1d1 )/(1j*omega*rho_e);
        vx_e1d2 =-derive_polynom_2D_x_2(p_e1d2 )/(1j*omega*rho_e);
        vx_e1d3 =-derive_polynom_2D_x_2(p_e1d3 )/(1j*omega*rho_e);
        vx_e1d4 =-derive_polynom_2D_x_2(p_e1d4 )/(1j*omega*rho_e);
        vx_e1d5 =-derive_polynom_2D_x_2(p_e1d5 )/(1j*omega*rho_e);
        vx_e1d6 =-derive_polynom_2D_x_2(p_e1d6 )/(1j*omega*rho_e);
        vx_e1d7 =-derive_polynom_2D_x_2(p_e1d7 )/(1j*omega*rho_e);
        vx_e1d8 =-derive_polynom_2D_x_2(p_e1d8 )/(1j*omega*rho_e);
        vx_e1d9 =-derive_polynom_2D_x_2(p_e1d9 )/(1j*omega*rho_e);
        vx_e1d10=-derive_polynom_2D_x_2(p_e1d10)/(1j*omega*rho_e);
        vx_e1d11=-derive_polynom_2D_x_2(p_e1d11)/(1j*omega*rho_e);
        vx_e1d12=-derive_polynom_2D_x_2(p_e1d12)/(1j*omega*rho_e);
        
        vy_e1d1 =-derive_polynom_2D_y_2(p_e1d1 )/(1j*omega*rho_e);
        vy_e1d2 =-derive_polynom_2D_y_2(p_e1d2 )/(1j*omega*rho_e);
        vy_e1d3 =-derive_polynom_2D_y_2(p_e1d3 )/(1j*omega*rho_e);
        vy_e1d4 =-derive_polynom_2D_y_2(p_e1d4 )/(1j*omega*rho_e);
        vy_e1d5 =-derive_polynom_2D_y_2(p_e1d5 )/(1j*omega*rho_e);
        vy_e1d6 =-derive_polynom_2D_y_2(p_e1d6 )/(1j*omega*rho_e);
        vy_e1d7 =-derive_polynom_2D_y_2(p_e1d7 )/(1j*omega*rho_e);
        vy_e1d8 =-derive_polynom_2D_y_2(p_e1d8 )/(1j*omega*rho_e);
        vy_e1d9 =-derive_polynom_2D_y_2(p_e1d9 )/(1j*omega*rho_e);
        vy_e1d10=-derive_polynom_2D_y_2(p_e1d10)/(1j*omega*rho_e);
        vy_e1d11=-derive_polynom_2D_y_2(p_e1d11)/(1j*omega*rho_e);
        vy_e1d12=-derive_polynom_2D_y_2(p_e1d12)/(1j*omega*rho_e);
        
        
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
        vx_e2d1=-derive_polynom_2D_x_2(p_e2d1)/(1j*omega*rho_e);
        vx_e2d2=-derive_polynom_2D_x_2(p_e2d2)/(1j*omega*rho_e);
        vx_e2d3=-derive_polynom_2D_x_2(p_e2d3)/(1j*omega*rho_e);
        vx_e2d4=-derive_polynom_2D_x_2(p_e2d4)/(1j*omega*rho_e);
        vx_e2d5=-derive_polynom_2D_x_2(p_e2d5)/(1j*omega*rho_e);
        vx_e2d6=-derive_polynom_2D_x_2(p_e2d6)/(1j*omega*rho_e);
        vy_e2d1=-derive_polynom_2D_y_2(p_e2d1)/(1j*omega*rho_e);
        vy_e2d2=-derive_polynom_2D_y_2(p_e2d2)/(1j*omega*rho_e);
        vy_e2d3=-derive_polynom_2D_y_2(p_e2d3)/(1j*omega*rho_e);
        vy_e2d4=-derive_polynom_2D_y_2(p_e2d4)/(1j*omega*rho_e);
        vy_e2d5=-derive_polynom_2D_y_2(p_e2d5)/(1j*omega*rho_e);
        vy_e2d6=-derive_polynom_2D_y_2(p_e2d6)/(1j*omega*rho_e);
    case 3
        index_p_2=dof_A(p_TR(elem.nodes(e_2,1:3)));
        nb_dof_2=3;
        vcor=nodes(nonzeros(elem.nodes(e_2,1:3)),1:2);
        p_e2d1=Lagrange_TR3(vcor,1);
        p_e2d2=Lagrange_TR3(vcor,2);
        p_e2d3=Lagrange_TR3(vcor,3);
        vx_e2d1=-derive_polynom_2D_x_2(p_e2d1)/(1j*omega*rho_e);
        vx_e2d2=-derive_polynom_2D_x_2(p_e2d2)/(1j*omega*rho_e);
        vx_e2d3=-derive_polynom_2D_x_2(p_e2d3)/(1j*omega*rho_e);
        vy_e2d1=-derive_polynom_2D_y_2(p_e2d1)/(1j*omega*rho_e);
        vy_e2d2=-derive_polynom_2D_y_2(p_e2d2)/(1j*omega*rho_e);
        vy_e2d3=-derive_polynom_2D_y_2(p_e2d3)/(1j*omega*rho_e);
    case 2
        index_p_2=dof_A(p_H12(elem.nodes(e_2,1:4)));
        nb_dof_2=12;
        
        lx=norm(nodes(elem.nodes(e_2,1),:)-nodes(elem.nodes(e_2,2),:));
        ly=norm(nodes(elem.nodes(e_2,1),:)-nodes(elem.nodes(e_2,4),:));
        [p_e2d1,p_e2d2,p_e2d3,p_e2d4,p_e2d5,p_e2d6,p_e2d7,p_e2d8,p_e2d9,p_e2d10,p_e2d11,p_e2d12]=H12_shape_functions_shifted(lx,ly,nodes(elem.nodes(e_2,1),1),nodes(elem.nodes(e_2,1),2));
        
        vx_e2d1 =-derive_polynom_2D_x_2(p_e2d1 )/(1j*omega*rho_e);
        vx_e2d2 =-derive_polynom_2D_x_2(p_e2d2 )/(1j*omega*rho_e);
        vx_e2d3 =-derive_polynom_2D_x_2(p_e2d3 )/(1j*omega*rho_e);
        vx_e2d4 =-derive_polynom_2D_x_2(p_e2d4 )/(1j*omega*rho_e);
        vx_e2d5 =-derive_polynom_2D_x_2(p_e2d5 )/(1j*omega*rho_e);
        vx_e2d6 =-derive_polynom_2D_x_2(p_e2d6 )/(1j*omega*rho_e);
        vx_e2d7 =-derive_polynom_2D_x_2(p_e2d7 )/(1j*omega*rho_e);
        vx_e2d8 =-derive_polynom_2D_x_2(p_e2d8 )/(1j*omega*rho_e);
        vx_e2d9 =-derive_polynom_2D_x_2(p_e2d9 )/(1j*omega*rho_e);
        vx_e2d10=-derive_polynom_2D_x_2(p_e2d10)/(1j*omega*rho_e);
        vx_e2d11=-derive_polynom_2D_x_2(p_e2d11)/(1j*omega*rho_e);
        vx_e2d12=-derive_polynom_2D_x_2(p_e2d12)/(1j*omega*rho_e);
        
        vy_e2d1 =-derive_polynom_2D_y_2(p_e2d1 )/(1j*omega*rho_e);
        vy_e2d2 =-derive_polynom_2D_y_2(p_e2d2 )/(1j*omega*rho_e);
        vy_e2d3 =-derive_polynom_2D_y_2(p_e2d3 )/(1j*omega*rho_e);
        vy_e2d4 =-derive_polynom_2D_y_2(p_e2d4 )/(1j*omega*rho_e);
        vy_e2d5 =-derive_polynom_2D_y_2(p_e2d5 )/(1j*omega*rho_e);
        vy_e2d6 =-derive_polynom_2D_y_2(p_e2d6 )/(1j*omega*rho_e);
        vy_e2d7 =-derive_polynom_2D_y_2(p_e2d7 )/(1j*omega*rho_e);
        vy_e2d8 =-derive_polynom_2D_y_2(p_e2d8 )/(1j*omega*rho_e);
        vy_e2d9 =-derive_polynom_2D_y_2(p_e2d9 )/(1j*omega*rho_e);
        vy_e2d10=-derive_polynom_2D_y_2(p_e2d10)/(1j*omega*rho_e);
        vy_e2d11=-derive_polynom_2D_y_2(p_e2d11)/(1j*omega*rho_e);
        vy_e2d12=-derive_polynom_2D_y_2(p_e2d12)/(1j*omega*rho_e);
        
        
end



%normal_displacement_e1=(v_xn_x+x_yn_y)/(j*omega)
temp=W_e1_in*Omega_e1_in;
Boundary_11= (temp(1,:)*nx+temp(2,:)*ny)/(1j*omega);
temp=W_e1_out*Omega_e2_in;
Boundary_12= (temp(1,:)*nx+temp(2,:)*ny)/(1j*omega);
%normal_displacement_e2=-(v_xn_x+x_yn_y)/(j*omega)
temp=W_e2_out*Omega_e1_in;
Boundary_21=-(temp(1,:)*nx+temp(2,:)*ny)/(1j*omega);
temp=W_e2_in*Omega_e2_in;
Boundary_22=-(temp(1,:)*nx+temp(2,:)*ny)/(1j*omega);



for i_test=1:nb_dof_1
    eval(['Interp_test=p_e1d',num2str(i_test),';']);
    for i_champs=1:nb_dof_1
        eval(['Interp_champs=Boundary_11(1)*vx_e1d',num2str(i_champs),'+Boundary_11(2)*vy_e1d',num2str(i_champs),'+Boundary_11(3)*p_e1d',num2str(i_champs),';'])
        A(index_p_1(i_test),index_p_1(i_champs))=A(index_p_1(i_test),index_p_1(i_champs))-integrate_polynom_2D_edge(multiply_polynom_2D(Interp_test,Interp_champs),a,b);
    end
    for i_champs=1:nb_dof_2
        eval(['Interp_champs=Boundary_12(1)*vx_e2d',num2str(i_champs),'+Boundary_12(2)*vy_e2d',num2str(i_champs),'+Boundary_12(3)*p_e2d',num2str(i_champs),';'])
        A(index_p_1(i_test),index_p_2(i_champs))=A(index_p_1(i_test),index_p_2(i_champs))-integrate_polynom_2D_edge(multiply_polynom_2D(Interp_test,Interp_champs),a,b);
    end
end

for i_test=1:nb_dof_2
    eval(['Interp_test=p_e2d',num2str(i_test),';']);
    for i_champs=1:nb_dof_1
        eval(['Interp_champs=Boundary_21(1)*vx_e1d',num2str(i_champs),'+Boundary_21(2)*vy_e1d',num2str(i_champs),'+Boundary_21(3)*p_e1d',num2str(i_champs),';'])
        A(index_p_2(i_test),index_p_1(i_champs))=A(index_p_2(i_test),index_p_1(i_champs))-integrate_polynom_2D_edge(multiply_polynom_2D(Interp_test,Interp_champs),a,b);
    end
    for i_champs=1:nb_dof_2
        eval(['Interp_champs=Boundary_22(1)*vx_e2d',num2str(i_champs),'+Boundary_22(2)*vy_e2d',num2str(i_champs),'+Boundary_22(3)*p_e2d',num2str(i_champs),';'])
        A(index_p_2(i_test),index_p_2(i_champs))=A(index_p_2(i_test),index_p_2(i_champs))-integrate_polynom_2D_edge(multiply_polynom_2D(Interp_test,Interp_champs),a,b);
    end
end

