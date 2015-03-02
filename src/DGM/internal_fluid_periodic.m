% internal_fluid_periodic.m
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
e_edge=e_left;
parameter_element

[F_plus,F_moins]=Split_fluid(nx,ny,Z_e);


% 
% B_e=[nx ny 0; 0 0 1];
% B_e_p=B_e/delta;
% 
% W_e_plus= Phi_fluid( nx, ny,Z_e);
% W_e_moins=Phi_fluid(-nx,-ny,Z_e);
% W_e_0=Phi_fluid_0(nx,ny);
% 
% 
% Omega_e=inv([W_e_plus W_e_moins W_e_0]);
% Omega_e_plus=Omega_e(1,:);
% Omega_e_moins=Omega_e(2,:);
% 
% Lambda_e_moins=omega/k_e;
% Lambda_e_plus=-omega/k_e;
% 
% W_e_p_plus= Phi_fluid(-nx,-ny,Z_e);
% W_e_p_moins=Phi_fluid( nx, ny,Z_e);
% W_e_p_0=Phi_fluid_0(nx,ny);
% 
% Omega_e_p=inv([W_e_p_plus W_e_p_moins W_e_p_0]);
% Omega_e_p_plus=Omega_e_p(1,:);
% Omega_e_p_moins=Omega_e_p(2,:);
% 
% 
% Lambda_e_p_moins=-omega/k_e;
% Lambda_e_p_plus=  omega/k_e;
% 
% M1=[-B_e*M_e*W_e_plus*Lambda_e_plus B_e_p*M_e*W_e_p_plus*Lambda_e_p_plus];
% M2=[ B_e*M_e*W_e_moins*Lambda_e_moins -B_e_p*M_e*W_e_p_moins*Lambda_e_p_moins];
% 
% Refl=inv(M2)*M1;
% 
% R_ee=Refl(1,1);
% R_eep=Refl(1,2);
% R_epe=Refl(2,1);
% R_epep=Refl(2,2);
% 
% 
% F_ee=  M_e*(W_e_plus*Lambda_e_plus+W_e_moins*Lambda_e_moins*R_ee)*Omega_e_plus;
% F_eep= M_e*W_e_moins*Lambda_e_moins*R_eep*Omega_e_p_plus;
% F_epe=-M_e*W_e_p_moins*Lambda_e_p_moins*R_epe*Omega_e_plus;
% F_epep=-M_e*(W_e_p_moins*Lambda_e_p_moins*R_epep+W_e_p_plus*Lambda_e_p_plus)*Omega_e_p_plus;
% 
% 
% for i_theta_test=1:nb_theta
%     theta_test=vec_theta(i_theta_test);
%     n_psi=[cos(theta_test);sin(theta_test)];
%     for  i_theta_champs=(1:nb_theta)
%         theta_champs=vec_theta(i_theta_champs);
%         n_phi=[cos(theta_champs);sin(theta_champs)];
%         % e'=e_1
%         % e=e_2
%         
%         % Champs direct et test
%         Psi_e=conj(Phi_fluid(cos(theta_test),sin(theta_test),Z_e)); 
%         Phi_e=     Phi_fluid(cos(theta_champs),sin(theta_champs),Z_e);
%         Psi_e_p=conj(Phi_fluid(cos(theta_test),sin(theta_test),Z_e)); 
%         Phi_e_p=     Phi_fluid(cos(theta_champs),sin(theta_champs),Z_e);
% 
%         
%         % ve^H F11 ue 
%         indice_test  =indice_fluid(e_right,i_theta_test,dof_start_element);
%         indice_champs=indice_fluid(e_right,i_theta_champs,dof_start_element);
%         
%         A(indice_test,indice_champs)=A(indice_test,indice_champs)+Psi_e'*F_ee *Phi_e...
%             *int_edge_2(j*k_e*n_psi,-j*k_e*n_phi,a_left,b_left,[c_right-[period;0] c_right-[period;0]]);
%         
%         % ve^H Feep ue'
%         
%         indice_test  =indice_fluid(e_left,i_theta_test,dof_start_element);
%         indice_champs=indice_fluid(e_right,i_theta_champs,dof_start_element);
% 
%         A(indice_test,indice_champs)=A(indice_test,indice_champs)+Psi_e'*F_eep *Phi_e_p...
%             *int_edge_2(j*k_e*n_psi,-j*k_e*n_phi,a_left,b_left,[c_left c_right-[period;0]]);
%         % -ve'^H F+ ue
%         
%         indice_test  =indice_fluid(e_right,i_theta_test,dof_start_element);
%         indice_champs=indice_fluid(e_left,i_theta_champs,dof_start_element);
%         
%        % F_epe
%         
%         
%         A(indice_test,indice_champs)=A(indice_test,indice_champs)+Psi_e_p'*F_epe *Phi_e...
%             *int_edge_2(j*k_e*n_psi,-j*k_e*n_phi,a_left,b_left,[c_right-[period;0] c_left]);
%         % -ve'^H F- ue'
%         
%         indice_test  =indice_fluid(e_left,i_theta_test,dof_start_element);
%         indice_champs=indice_fluid(e_left,i_theta_champs,dof_start_element);
%         A(indice_test,indice_champs)=A(indice_test,indice_champs)+Psi_e_p'*F_epep *Phi_e_p...
%             *int_edge_2(j*k_e*n_psi,-j*k_e*n_phi,a_left,b_left,[c_left c_left]);
%     end
% end



%%All the integrations are on the left boundary

nx=cos(vec_theta);
ny=sin(vec_theta);
Phi=Phi_fluid_vector(nx,ny,Z_e,Shift_fluid);

II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a_left,b_left,[c_right-[period;0] c_right-[period;0]]);
MM=kron(II,F_plus);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_right);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_right);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi/delta;

        
II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a_left,b_left,[c_right-[period;0] c_left]);
MM=kron(II,F_moins);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_right);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_left);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi/delta;

II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a_left,b_left,[c_left c_right-[period;0]]);
MM=kron(II,-F_plus);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_left);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_right);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;

II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a_left,b_left,[c_left c_left]);
MM=kron(II,-F_moins);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_left);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_left);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;

