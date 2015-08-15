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

e_1=edges.internal(ie,3);
e_2=edges.internal(ie,4);
c_1=mean(nodes(nonzeros(elem.nodes(e_1,:)),1:2))';
c_2=mean(nodes(nonzeros(elem.nodes(e_2,:)),1:2))';
e_edge=e_1;


coord_edge(1:2,1)=nodes(edges.internal(ie,1),1:2)';
coord_edge(1:2,2)=nodes(edges.internal(ie,2),1:2)';

a=coord_edge(:,1);
b=coord_edge(:,2);

h=norm(b-a);
n_ab=(b-a)/h;

%%%%% Elements on both sides of the edge



%%%%% vector normal pointing out from e_1

centre_edge=(a+b)/2;
n_centre=c_1-centre_edge;
ne=normal_edge(coord_edge);
if (n_centre'*ne>0)
    ne=-ne;
end
nx=ne(1);
ny=ne(2);


% % u=[vx vy p]
% %i*omega*u+A_x du/dx+A_y dudy
%
A_x=[0 1/air.rho 0;0 0 0; air.K 0 0];
A_y=[0 0 0;0 0 1/air.rho; 0 air.K 0];
F_e=-M_e*(A_x*nx+A_y*ny);
% %


W_e1_in=[nx;ny;air.Z];
Omega_e1_in=[nx/2 ny/2 1/(2*air.Z)];
W_e1_out=[nx;ny;-air.Z];
Omega_e1_out=[nx/2 ny/2 -1/(2*air.Z)];

W_e2_in=W_e1_out;
Omega_e2_in=Omega_e1_out;
W_e2_out=W_e1_in;
Omega_e2_out=Omega_e1_in;


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

indice_DGM  =((1:theta_DGM.nb)-1)+dof_start_element(e_2);
 
%normal_displacement_e1=(v_xn_x+x_yn_y)/(j*omega)
temp=W_e1_in*Omega_e1_in;
Boundary_11= (temp(1,:)*nx+temp(2,:)*ny)/(1j*omega);
temp=W_e1_out*Omega_e2_in;
Boundary_12= (temp(1,:)*nx+temp(2,:)*ny)/(1j*omega);


for i_test=1:nb_dof_1
    eval(['Interp_test=p_e1d',num2str(i_test),';']);
    for i_champs=1:nb_dof_1
        eval(['Interp_champs=Boundary_11(1)*vx_e1d',num2str(i_champs),'+Boundary_11(2)*vy_e1d',num2str(i_champs),'+Boundary_11(3)*p_e1d',num2str(i_champs),';'])
        A(index_p_1(i_test),index_p_1(i_champs))=A(index_p_1(i_test),index_p_1(i_champs))-integrate_polynom_2D_edge(multiply_polynom_2D(Interp_test,Interp_champs),a,b);
    end
    for i_champs=1:theta_DGM.nb
        delta_exp=int_edge_1vectorielle(-1j*k_e*[cos(vec_theta(i_champs));sin(vec_theta(i_champs))],a,b,c_2)/norm(b-a);
        Phi_champs=Phi_fluid(cos(vec_theta(i_champs)),sin(vec_theta(i_champs)),Z_e);
        temp=Boundary_12*Phi_champs;
        A(index_p_1(i_test),indice_DGM(i_champs))=A(index_p_1(i_test),indice_DGM(i_champs))-delta_exp*temp*integrate_polynom_2D_edge(Interp_test,a,b);
    end
end

for i_test=1:theta_DGM.nb
    Phi_test=Phi_fluid(cos(vec_theta(i_test)),sin(vec_theta(i_test)),Z_e);
    for i_champs=1:nb_dof_1
        delta_exp=int_edge_1vectorielle(1j*k_e*[cos(vec_theta(i_test));sin(vec_theta(i_test))],a,b,c_2)/norm(b-a);
        temp=Phi_test.'*W_e2_out*Omega_e1_in;
        eval(['Interp_champs=temp(1)*vx_e1d',num2str(i_champs),'+temp(2)*vy_e1d',num2str(i_champs),'+temp(3)*p_e1d',num2str(i_champs),';'])
        A(indice_DGM(i_test),index_p_1(i_champs))=A(indice_DGM(i_test),index_p_1(i_champs))-delta_exp*integrate_polynom_2D_edge(Interp_champs,a,b);
    end
end




FF=F_e*W_e2_in*Omega_e2_in;
nx=cos(vec_theta);
ny=sin(vec_theta);
Phi=Phi_fluid_vector(nx,ny,Z_e,Shift_fluid);
II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a,b,[c_2 c_2]);
MM=kron(II,FF);
A(indice_DGM,indice_DGM)=A(indice_DGM,indice_DGM)-Phi.'*MM*Phi;
