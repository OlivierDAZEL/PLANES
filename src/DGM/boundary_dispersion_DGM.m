% boundary_dispersion_DGM.m
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


n_excitation=[cos(data_model.theta_inc); sin(data_model.theta_inc)];

e_edge=edges.loads(ie,3);
c_e= mean(nodes(nonzeros(elem.nodes(e_edge,:)),1:2))';
coord_edge(1:2,1)=nodes(edges.loads(ie,1),1:2)';
coord_edge(1:2,2)=nodes(edges.loads(ie,2),1:2)';
a=coord_edge(:,1);
b=coord_edge(:,2);

[nx, ny] = normal_edge_out_element(a, b, c_e);

if (elem.label(e_edge)==0)

    M_e=diag([air.rho,air.rho,1/air.K]);

    A_x=[0 0 1/air.rho;0 0 0; air.rho*air.c^2 0 0];
    A_y=[0 0 0;0 0 1/air.rho; 0 air.rho*air.c^2 0];

    F_e=(A_x*nx+A_y*ny);
		[P_e, Lambda_e]  = eig(F_e);
		Q_e = inv(P_e);

		% upwinding
		rlp_list = find(diag(Lambda_e)<=0);
		rlm_list = find(diag(Lambda_e)>=0);

		lambda_minus = Lambda_e; lambda_minus(:,rlm_list) = 0;
		lambda_plus = Lambda_e; lambda_plus(:,rlp_list) = 0;

		P_e_in = P_e; P_e_in(:,rlp_list) = 0;
		P_e_out = P_e; P_e_out(:,rlm_list) = 0;

		Q_e_in = Q_e; Q_e_in(rlp_list,:) = 0;
		Q_e_out = Q_e; Q_e_out(rlm_list,:) = 0;

		F_e_in = P_e*lambda_plus*Q_e;
		F_e_out = P_e*lambda_minus*Q_e;

    nx=cos(vec_theta);
    ny=sin(vec_theta);
    indices  =((1:data_model.theta_DGM.nb)-1)+dof_start_element(e_edge);
    Phi=Phi_fluid_vector(nx, ny, air.Z, Shift_fluid);
    k_e = omega/air.c;

    % Modified Flux Matrix
    II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a,b,[c_e c_e]);
    A(indices,indices)=A(indices,indices)+Phi.'*kron(II,M_e*F_e_in)*Phi;

    % Forcing term
    % II=int_edge_1vectorielle_test(1j*k_e*[nx;ny],-1j*k_e*[n_excitation(1)*ones(1,length(nx)); n_excitation(2)*ones(1,length(ny))],a,b,c_e);
    II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*n_excitation,a,b,[c_e, [0;0]]);
    Phi_excit=Phi_fluid(n_excitation(1)/air.Z, n_excitation(2)/air.Z, 1);
    F(indices)=F(indices)-Phi.'*kron(II, M_e*F_e_out)*Phi_excit;


end
