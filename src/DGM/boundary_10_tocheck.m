% boundary_10.m
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
%%%%% Coordinates of the edge's vertices

n_excitation=[sin(theta_inc); cos(theta_inc)];

coord_edge(1:2,1)=nodes(loads(ie,1),1:2)';
coord_edge(1:2,2)=nodes(loads(ie,2),1:2)';

a=coord_edge(:,1);
b=coord_edge(:,2);

h=norm(b-a);
n_ab=(b-a)/h;

%%%%% Element linked to the edge

e_edge=loads(ie,3);
c_edge=mean(nodes(elements(e_edge,:),1:2))';

parameter_element

%%%%% vector normal pointing towards \Omega_e

centre_temp=mean(nodes(elements(e_edge,:),1:2))';
centre_edge=(a+b)/2;
n_centre=centre_temp-centre_edge;
ne=normal_edge(coord_edge);
if (n_centre'*ne>0)
    ne=-ne;
end
nx=ne(1);
ny=ne(2);

if (element_label(e_edge)==0)
    
    W_e_plus=Phi_fluid(nx,ny,Z_e);
    W_e_moins=Phi_fluid(-nx,-ny,Z_e);
    W_e_0=Phi_fluid_0(nx,ny);
    
    
    W_e=[W_e_plus W_e_moins W_e_0];
    Omega_e=inv(W_e);
    
    Omega_e_plus=Omega_e(1,:);
    Omega_e_moins=Omega_e(2,:);
    
    
    Lambda_e_plus= diag(c_e);
    Lambda_e_moins=-Lambda_e_plus;
    
    B_e=[0 0 1];
    s_e=-cos(theta_inc)/Z_e;
    
    M_plus =B_e*M_e*W_e_plus *Lambda_e_plus;
    M_moins=B_e*M_e*W_e_moins*Lambda_e_moins;
    
    
    R_e=-inv(M_moins)*M_plus;
    S_e= inv(M_moins)*s_e;
    
    
    %F_e=A_x_fluid*nx+A_y_fluid*ny;
    F_e=M_e*(Lambda_e_plus*W_e_plus+Lambda_e_moins*W_e_moins*R_e)*Omega_e_plus;
    S_e=M_e*                        Lambda_e_moins*W_e_moins*S_e;
    
    
    for i_thetapsi=1:nb_thetaDGM
        theta_psi=vec_theta(i_thetapsi);
        n_psi=[cos(theta_psi);sin(theta_psi)];
        Psi_e=conj(Phi_fluid(cos(theta_psi),sin(theta_psi),Z_e));
        ii=indice_fluid(e_edge,i_thetapsi,dof_start_element);
        
        for i_thetaphi=(1:nb_thetaDGM)
            theta_phi=vec_theta(i_thetaphi);
            n_phi=[cos(theta_phi);sin(theta_phi)];
            Phi_e=     Phi_fluid(cos(theta_phi),sin(theta_phi),Z_e) ;
            jj=indice_fluid(e_edge,i_thetaphi,dof_start_element);
            A(ii,jj)=A(ii,jj)+Psi_e'*F_e*Phi_e*...
                int_edge_2(1j*k_e*n_psi,-1j*k_e*n_phi,a,b,[c_edge c_edge]);
        end
    end
    
    
    for i_thetapsi=1:nb_thetaDGM
        theta_psi=vec_theta(i_thetapsi);
        n_psi=[cos(theta_psi);sin(theta_psi)];
        Psi_e=conj(Phi_fluid(cos(theta_psi),sin(theta_psi),Z_e));
        ii=indice_fluid(e_edge,i_thetapsi,dof_start_element);
        F(ii)=F(ii)-Psi_e'*S_e*...
            int_edge_2(1j*k_e*n_psi,-1j*k_e*n_excitation,a,b,[c_edge 0*c_edge]);
        
        A(ii,nb.dof_DGM+1)=A(ii,nb.dof_DGM+1)-Psi_e'*S_e*...
            int_edge_2(1j*k_e*n_psi,-1j*k_e*n_excitation,a,b,[c_edge 0*c_edge]);
    end
    
    
    for i_thetaphi=(1:nb_thetaDGM)
        theta_phi=vec_theta(i_thetaphi);
        n_phi=[cos(theta_phi);sin(theta_phi)];
        Phi_e=     Phi_fluid(cos(theta_phi),sin(theta_phi),Z_e) ;
        jj=indice_fluid(e_edge,i_thetaphi,dof_start_element);
        A(nb.dof_DGM+1,jj)=A(nb.dof_DGM+1,jj)+Phi_e(3)*...
            int_edge_2(1j*k_e*n_excitation,-1j*k_e*n_phi,a,b,[0*c_edge c_edge]);
    end
    A(nb.dof_DGM+1,nb.dof_DGM+1)=-period;
    F(nb.dof_DGM+1)=period;
    
end


if (floor(element_label(e_edge)/1000)==4)
    
    W_e_plus= Phi_Biot( nx, ny,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
    W_e_moins=Phi_Biot(-nx,-ny,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
    W_e_0=Phi_Biot_0(nx,ny);
    
    
    W_e=[W_e_plus W_e_moins W_e_0];
    Omega_e=inv(W_e);
    
    Omega_e_plus=Omega_e(1:3,:);
    Omega_e_moins=Omega_e(4:6,:);
    
    Lambda_e_plus=diag([omega/delta_1 omega/delta_2 omega/delta_3]);
    Lambda_e_moins=-Lambda_e_plus;
    
    B_e=zeros(3,8);
    s_e=zeros(3,1);
    
    B_e(1,1)=1;
    B_e(2,2)=1;
    B_e(3,8)=1;
    s_e(3)=-cos(theta_inc)/air.Z;
    
    M_plus =B_e*M_e*W_e_plus *Lambda_e_plus;
    M_moins=B_e*M_e*W_e_moins*Lambda_e_moins;
    
    
    R_e=-inv(M_moins)*M_plus;
    S_e= inv(M_moins)*s_e;
    
    
    %F_e=A_x_fluid*nx+A_y_fluid*ny;
    F_e=M_e*(W_e_plus*Lambda_e_plus+W_e_moins*Lambda_e_moins*R_e)*Omega_e_plus;
    S_e=M_e*                        W_e_moins*Lambda_e_moins*S_e;
    
    
    nx=cos(vec_theta);
    ny=sin(vec_theta);
    Phi=Phi_Biot_vector(nx,ny,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega,Shift_Biot);
    
    II=int_edge_2vectorielle(1j*[delta_1*[nx;ny] delta_2*[nx;ny] delta_3*[nx;ny]],-1j*[delta_1*[nx;ny] delta_2*[nx;ny] delta_3*[nx;ny]],a,b,[c_edge c_edge]);
    MM=kron(II,F_e);
    indice_test  =[1+3*(0:nb_thetaDGM-1) 2+3*(0:nb_thetaDGM-1) 3+3*(0:nb_thetaDGM-1)]+dof_start_element(e_edge)-1;
    indice_champs=[1+3*(0:nb_thetaDGM-1) 2+3*(0:nb_thetaDGM-1) 3+3*(0:nb_thetaDGM-1)]+dof_start_element(e_edge)-1;
    A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;
    
    
    II=int_edge_2vectorielle(1j*[delta_1*[nx;ny] delta_2*[nx;ny] delta_3*[nx;ny]],-1j*k_air*n_excitation,a,b,[c_edge 0*c_edge]);
    MM=kron(II,S_e);
    indice_test  =[1+3*(0:nb_thetaDGM-1) 2+3*(0:nb_thetaDGM-1) 3+3*(0:nb_thetaDGM-1)]+dof_start_element(e_edge)-1;
    F(indice_test)=F(indice_test)-Phi.'*MM;
    A(indice_test,nb.dof_DGM+1)=A(indice_test,nb.dof_DGM+1)-Phi.'*MM;
    
    II=int_edge_2vectorielle(1j*k_air*n_excitation,-1j*[delta_1*[nx;ny] delta_2*[nx;ny] delta_3*[nx;ny]],a,b,[0*c_edge c_edge]);
    
    % Matrix to filter sigma_p
    S_filter=zeros(1,3*nb_thetaDGM*8);
    S_filter(1,8*[1:3*nb_thetaDGM])=1;
    Phi=S_filter*Phi;
    
    indice_champs  =[1+3*(0:nb_thetaDGM-1) 2+3*(0:nb_thetaDGM-1) 3+3*(0:nb_thetaDGM-1)]+dof_start_element(e_edge)-1;
    A(nb.dof_DGM+1,indice_champs)=A(nb.dof_DGM+1,indice_champs)+Phi.*II;

    A(nb.dof_DGM+1,nb.dof_DGM+1)=-period;
    F(nb.dof_DGM+1)=period;
    
end



