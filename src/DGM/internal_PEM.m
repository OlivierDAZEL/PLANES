% internal_PEM.m
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

parameter_element

W_e_plus= Phi_Biot(-nx,-ny,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
W_e_moins=Phi_Biot( nx, ny,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
W_e_0=Phi_Biot_0(nx,ny);

Omega_e=inv([W_e_plus W_e_moins W_e_0]);

Omega_e_plus=Omega_e(1:3,:);
Omega_e_moins=Omega_e(4:6,:);

Lambda_e_moins=-diag([omega/delta_1 omega/delta_2 omega/delta_3]);
Lambda_e_plus=-Lambda_e_moins;

F_plus=M_e*W_e_plus*Lambda_e_plus*Omega_e_plus;
F_moins=M_e*W_e_moins*Lambda_e_moins*Omega_e_moins;

%F_plus+F_moins



delta=[delta_1 delta_2 delta_3];


for i_thetapsi=1:nb_theta
    theta_test=vec_theta(i_thetapsi);
    n_psi=[cos(theta_test);sin(theta_test)];
    Psi_e=conj(Phi_Biot(cos(theta_test)  ,sin(theta_test)  ,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega));
    for  i_thetaphi=(1:nb_theta)
        theta_champs=vec_theta(i_thetaphi);
        n_phi=[cos(theta_champs);sin(theta_champs)];

        Phi_e=     Phi_Biot(cos(theta_champs),sin(theta_champs),delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
        i_test=1:3; % Balayage des ondes de Biot champs test
        i_champs=1:3;

        
        integrale_border11=int_edge_2vectorielle(j*(ones(2,1)*delta).*(n_psi*ones(1,3)),-j*(ones(2,1)*delta).*(n_phi*ones(1,3)),a,b,[c_1 c_1]);
        integrale_border12=int_edge_2vectorielle(j*(ones(2,1)*delta).*(n_psi*ones(1,3)),-j*(ones(2,1)*delta).*(n_phi*ones(1,3)),a,b,[c_1 c_2]);
        integrale_border21=int_edge_2vectorielle(j*(ones(2,1)*delta).*(n_psi*ones(1,3)),-j*(ones(2,1)*delta).*(n_phi*ones(1,3)),a,b,[c_2 c_1]);
        integrale_border22=int_edge_2vectorielle(j*(ones(2,1)*delta).*(n_psi*ones(1,3)),-j*(ones(2,1)*delta).*(n_phi*ones(1,3)),a,b,[c_2 c_2]);

        
        %ve^H F+ ue
        ii=indice_Biot(e_2,i_thetapsi,i_test,dof_start_element);
        jj=indice_Biot(e_2,i_thetaphi,i_champs,dof_start_element);
        A(ii,jj)=A(ii,jj)+(Psi_e(:,i_test)'*F_plus *Phi_e(:,i_champs)).*integrale_border22;
        %ve^H F- ue'
        ii=indice_Biot(e_2,i_thetapsi,i_test,dof_start_element);
        jj=indice_Biot(e_1,i_thetaphi,i_champs,dof_start_element);
        A(ii,jj)=A(ii,jj)+Psi_e(:,i_test)'*F_moins *Phi_e(:,i_champs).*integrale_border21;
        %-ve'^H F+ ue
        ii=indice_Biot(e_1,i_thetapsi,i_test,dof_start_element);
        jj=indice_Biot(e_2,i_thetaphi,i_champs,dof_start_element);
        A(ii,jj)=A(ii,jj)-Psi_e(:,i_test)'*F_plus *Phi_e(:,i_champs).*integrale_border12;
        %-ve'^H F- ue'
        ii=indice_Biot(e_1,i_thetapsi,i_test,dof_start_element);
        jj=indice_Biot(e_1,i_thetaphi,i_champs,dof_start_element);
        A(ii,jj)=A(ii,jj)-Psi_e(:,i_test)'*F_moins *Phi_e(:,i_champs).*integrale_border11;




    end
end




% for i_thetapsi=1:nb_theta
%     theta_test=vec_theta(i_thetapsi);
%     n_psi=[cos(theta_test);sin(theta_test)];
%     for  i_thetaphi=(1:nb_theta)
%         theta_champs=vec_theta(i_thetaphi);
%         n_phi=[cos(theta_champs);sin(theta_champs)];
%         Psi_e=conj(Phi_Biot(cos(theta_test)  ,sin(theta_test)  ,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega));
%         Phi_e=     Phi_Biot(cos(theta_champs),sin(theta_champs),delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
%         
%         integrale_border22=int_edge_2vectorielle(j*(ones(2,1)*delta).*(n_psi*ones(1,3)),-j*(ones(2,1)*delta).*(n_phi*ones(1,3)),a,b,[c_2 c_2])
%         
%         
%         for i_test=1:3 % Balayage des ondes de Biot champs test
%             for i_champs=1:3 % Balayage des ondes de Biot champs inconnu
%                 % ve^H F+ ue
%                 ii=indice_Biot(e_2,i_thetapsi,i_test);
%                 jj=indice_Biot(e_2,i_thetaphi,i_champs);
%                 M_DGM(ii,jj)=M_DGM(ii,jj)+Psi_e(:,i_test)'*F_plus *Phi_e(:,i_champs)...
%                         *int_edge_2(j*delta(i_test)*n_psi,-j*delta(i_champs)*n_phi,a,b,[c_2 c_2]);
%                 % ve^H F- ue'
%                 ii=indice_Biot(e_2,i_thetapsi,i_test);
%                 jj=indice_Biot(e_1,i_thetaphi,i_champs);
%                 M_DGM(ii,jj)=M_DGM(ii,jj)+Psi_e(:,i_test)'*F_moins *Phi_e(:,i_champs)...
%                         *int_edge_2(j*delta(i_test)*n_psi,-j*delta(i_champs)*n_phi,a,b,[c_2 c_1]);
%                 % -ve'^H F+ ue
%                 ii=indice_Biot(e_1,i_thetapsi,i_test);
%                 jj=indice_Biot(e_2,i_thetaphi,i_champs);
%                 M_DGM(ii,jj)=M_DGM(ii,jj)-Psi_e(:,i_test)'*F_plus *Phi_e(:,i_champs)...
%                         *int_edge_2(j*delta(i_test)*n_psi,-j*delta(i_champs)*n_phi,a,b,[c_1 c_2]);
%                 % -ve'^H F- ue'
%                 ii=indice_Biot(e_1,i_thetapsi,i_test);
%                 jj=indice_Biot(e_1,i_thetaphi,i_champs);
%                 M_DGM(ii,jj)=M_DGM(ii,jj)-Psi_e(:,i_test)'*F_moins *Phi_e(:,i_champs)...
%                     *int_edge_2(j*delta(i_test)*n_psi,-j*delta(i_champs)*n_phi,a,b,[c_1 c_1]);
% 
%             end
%         end
%     end
% end