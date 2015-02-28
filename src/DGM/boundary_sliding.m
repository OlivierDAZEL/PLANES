% boundary_sliding.m
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
coord_edge(1:2,1)=nodes(dirichlets(ie,1),1:2)';
coord_edge(1:2,2)=nodes(dirichlets(ie,2),1:2)';

a=coord_edge(:,1);
b=coord_edge(:,2);

h=norm(b-a);
n_ab=(b-a)/h;

%%%%% Element linked to the edge

e_edge=dirichlets(ie,3);
c_edge=mean(nodes(elements(e_edge,:),1:2))';

parameter_element


%%%%% vector normal pointing towards \Omega_e

centre_temp=mean(nodes(elements(e_edge,:),1:2))';
centre_edge=(a+b)/2;
n_centre=centre_temp-centre_edge;
ne=normal_edge(coord_edge);
if (n_centre'*ne<0)
    ne=-ne;
end
nx=ne(1);
ny=ne(2);


W_e_plus= Phi_Biot( nx, ny,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
W_e_moins=Phi_Biot(-nx,-ny,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
W_e_0=Phi_Biot_0(nx,ny);


W_e=[W_e_plus W_e_moins W_e_0];
Omega_e=inv(W_e);

Omega_e_plus=Omega_e(1:3,:);
Omega_e_moins=Omega_e(4:6,:);



delta=[delta_1 delta_2 delta_3];

Lambda_e_plus=diag([omega/delta_1 omega/delta_2 omega/delta_3]);
Lambda_e_moins=-Lambda_e_plus;

B_e=[-ny nx 0 0 0 0 0 0;0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 1];

M_plus =B_e*M_e*W_e_plus *Lambda_e_plus;
M_moins=B_e*M_e*W_e_moins*Lambda_e_moins;

R_e=-inv(M_moins)*M_plus;

%F_e=A_x_fluid*nx+A_y_fluid*ny;
F_e=M_e*(W_e_plus*Lambda_e_plus+W_e_moins*Lambda_e_moins*R_e)*Omega_e_plus;


for i_thetapsi=1:nb_theta
    theta_test=vec_theta(i_thetapsi);
    n_psi=[cos(theta_test);sin(theta_test)];
    for  i_thetaphi=(1:nb_theta)
        theta_champs=vec_theta(i_thetaphi);
        n_phi=[cos(theta_champs);sin(theta_champs)];
        Psi_e=conj(Phi_Biot(cos(theta_test),sin(theta_test),delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega));
        Phi_e=Phi_Biot(cos(theta_champs),sin(theta_champs),delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
        for i_test=1:3
            for i_champs=1:3
                % ve^H F+ ue
                ii=indice_Biot(e_edge,i_thetapsi,i_test,dof_start_element);
                jj=indice_Biot(e_edge,i_thetaphi,i_champs,dof_start_element);
                A(ii,jj)=A(ii,jj)+Psi_e(:,i_test)'*F_e*Phi_e(:,i_champs)...
                    *int_edge_2(j*delta(i_test)*n_psi,-j*delta(i_champs)*n_phi,a,b,[c_edge c_edge]);
            end
        end
    end
end