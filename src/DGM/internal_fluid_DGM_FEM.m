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
e_edge=e_1;
parameter_element

c_edge=mean(nodes(elements(e_edge,:),1:2))';

nx=0;
ny=1;

F_e=(M_e)*[0 0 0;0 0 1/air.rho; 0 air.K 0]*ny


C_e2=[0 1 0;0 0 1];
C_e1=[0 1 0;0 0 1];
W_e2=[nx -ny -nx;ny nx -ny;air.Z 0 air.Z];
Omega_e2=[nx/2 ny/2 1/(2*air.Z);-ny nx 0;-nx/2 -ny/2 1/(2*air.Z)];
W_e2_plus=W_e2(:,1);
W_e2_moins=W_e2(:,3);
Omega_e2_plus=Omega_e2(1,:);
Omega_e2_moins=Omega_e2(3,:);

W_e1=W_e2;
Omega_e1=Omega_e2;

W_e1_moins=W_e2_plus;
W_e1_plus=W_e2_moins;
Omega_e1_moins=Omega_e2_plus;
Omega_e1_plus=Omega_e2_moins;

R_tilde=inv([C_e2*W_e2_moins -C_e1*W_e1_plus])*([C_e1*W_e1_moins -C_e2*W_e2_plus]);


load_Hermite_2D


index_p=dof_A(p_H12(elements(e_2,:)));
indice_FEM_0=index_p([1 2 5 4]);
indice_FEM_y=index_p([3 6]);

index_Psi=[1 4];

indice_DGM  =((1:nb_thetaDGM)-1)+dof_start_element(e_1);

PP=W_e2_plus*Omega_e2_plus;
line_PP=2;
PP=-PP(line_PP,:)/(1j*omega); %\d p/\p n =-v_y
for i_test=1:4
    eval(['Psi_test=Psi_',num2str(i_test),'_x;'])
    for i_champs=1:4
        eval(['temp_1=PP(3)*Psi_',num2str(i_champs),'_x;']);
        eval(['temp_2=PP(1)*(-1/(j*air.rho*omega))*derive_polynom(Psi_',num2str(i_champs),'_x);']);
        temp=add_polynom(temp_1,temp_2);
        A(indice_FEM_0(i_test),indice_FEM_0(i_champs))=A(indice_FEM_0(i_test),indice_FEM_0(i_champs))-integrate_polynom(multiply_polynom(Psi_test,temp),lx_H12);
    end
    for i_champs=1:2
        eval(['temp=(1/ly_H12)*PP(2)*(-1/(j*air.rho*omega))*Psi_',num2str(index_Psi(i_champs)),'_x;'])
        A(indice_FEM_0(i_test),indice_FEM_y(i_champs))=A(indice_FEM_0(i_test),indice_FEM_y(i_champs))-integrate_polynom(multiply_polynom(Psi_test,temp),lx_H12);
    end
end


PP=W_e2_moins*Omega_e1_moins;
line_PP=2;
PP=-PP(line_PP,:)/(1j*omega); %\d p/\p n =-v_y

for i_test=1:4
    eval(['Psi_test=Psi_',num2str(i_test),'_x;']);
    for i_champs=1:nb_thetaDGM
        delta_exp=int_edge_1vectorielle(-1j*k_e*[cos(vec_theta(i_champs));sin(vec_theta(i_champs))],a,b,c_edge)/norm(b-a);
        Phi_champs=Phi_fluid(cos(vec_theta(i_champs)),sin(vec_theta(i_champs)),Z_e);
        temp=PP*Phi_champs;
        A(indice_FEM_0(i_test),indice_DGM(i_champs))=A(indice_FEM_0(i_test),indice_DGM(i_champs))-delta_exp*temp*integrate_polynom(Psi_test,lx_H12);
    end
end

FF=F_e*W_e1_plus*Omega_e2_plus;
for i_test=1:nb_thetaDGM
    delta_exp=int_edge_1vectorielle(1j*k_e*[cos(vec_theta(i_test));sin(vec_theta(i_test))],a,b,c_edge)/norm(b-a);
    Phi_test=Phi_fluid(cos(vec_theta(i_test)),sin(vec_theta(i_test)),Z_e);
    temp=Phi_test.'*FF;
    for i_champs=1:4
        eval(['Psi_champs=(-1/(j*omega*air.rho))*derive_polynom(Psi_',num2str(i_champs),'_x);']);
        A(indice_DGM(i_test),indice_FEM_0(i_champs))=A(indice_DGM(i_test),indice_FEM_0(i_champs))+delta_exp*temp(1)*integrate_polynom(Psi_champs,lx_H12);
        eval(['Psi_champs=Psi_',num2str(i_champs),'_x;']);
        A(indice_DGM(i_test),indice_FEM_0(i_champs))=A(indice_DGM(i_test),indice_FEM_0(i_champs))+delta_exp*temp(3)*integrate_polynom(Psi_champs,lx_H12);
    end
    for i_champs=1:2
        eval(['Psi_champs=(1/ly_H12)*(-1/(j*air.rho*omega))*Psi_',num2str(index_Psi(i_champs)),'_x;'])
        A(indice_DGM(i_test),indice_FEM_y(i_champs))=A(indice_DGM(i_test),indice_FEM_y(i_champs))+delta_exp*temp(2)*integrate_polynom(Psi_champs,lx_H12);
    end
end


FF=F_e*W_e1_moins*Omega_e1_moins;
nx=cos(vec_theta);
ny=sin(vec_theta);
Phi=Phi_fluid_vector(nx,ny,Z_e,Shift_fluid);
II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a,b,[c_1 c_1]);
MM=kron(II,FF);
A(indice_DGM,indice_DGM)=A(indice_DGM,indice_DGM)+Phi.'*MM*Phi;
