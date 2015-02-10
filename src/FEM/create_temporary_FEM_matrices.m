% create_temporary_FEM_matrices.m
%
% Copyright (C) 2014 < Olivier DAZEL <olivier.dazel@univ-lemans.fr> >
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
% https://github.com/OlivierDAZEL/PLANES
% or find more details on Olivier's webpage
% http://perso.univ-lemans.fr/~odazel/
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
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%%


i_k0_pem01=zeros(1,nb_media.pem01);
j_k0_pem01=zeros(1,nb_media.pem01);
v_k0_pem01=zeros(1,nb_media.pem01);
i_k1_pem01=zeros(1,nb_media.pem01);
j_k1_pem01=zeros(1,nb_media.pem01);
v_k1_pem01=zeros(1,nb_media.pem01);
i_m_pem01=zeros(1,nb_media.pem01);
j_m_pem01=zeros(1,nb_media.pem01);
v_m_pem01=zeros(1,nb_media.pem01);
i_h_pem01=zeros(1,nb_media.pem01);
j_h_pem01=zeros(1,nb_media.pem01);
v_h_pem01=zeros(1,nb_media.pem01);
i_q_pem01=zeros(1,nb_media.pem01);
j_q_pem01=zeros(1,nb_media.pem01);
v_q_pem01=zeros(1,nb_media.pem01);
i_c_pem01=zeros(1,nb_media.pem01);
j_c_pem01=zeros(1,nb_media.pem01);
v_c_pem01=zeros(1,nb_media.pem01);
i_cp_pem01=zeros(1,nb_media.pem01);
j_cp_pem01=zeros(1,nb_media.pem01);
v_cp_pem01=zeros(1,nb_media.pem01);



i_k0_pem98=zeros(1,nb_media.pem98);
j_k0_pem98=zeros(1,nb_media.pem98);
v_k0_pem98=zeros(1,nb_media.pem98);
i_k1_pem98=zeros(1,nb_media.pem98);
j_k1_pem98=zeros(1,nb_media.pem98);
v_k1_pem98=zeros(1,nb_media.pem98);
i_m_pem98=zeros(1,nb_media.pem98);
j_m_pem98=zeros(1,nb_media.pem98);
v_m_pem98=zeros(1,nb_media.pem98);
i_h_pem98=zeros(1,nb_media.pem98);
j_h_pem98=zeros(1,nb_media.pem98);
v_h_pem98=zeros(1,nb_media.pem98);
i_q_pem98=zeros(1,nb_media.pem98);
j_q_pem98=zeros(1,nb_media.pem98);
v_q_pem98=zeros(1,nb_media.pem98);
i_c_pem98=zeros(1,nb_media.pem98);
j_c_pem98=zeros(1,nb_media.pem98);
v_c_pem98=zeros(1,nb_media.pem98);

i_k0_elas=zeros(1,nb_media.elas);
j_k0_elas=zeros(1,nb_media.elas);
v_k0_elas=zeros(1,nb_media.elas);
i_k1_elas=zeros(1,nb_media.elas);
j_k1_elas=zeros(1,nb_media.elas);
v_k1_elas=zeros(1,nb_media.elas);
i_m_elas=zeros(1,nb_media.elas);
j_m_elas=zeros(1,nb_media.elas);
v_m_elas=zeros(1,nb_media.elas);


i_h_acoustic=zeros(1,nb_media.acou);
j_h_acoustic=zeros(1,nb_media.acou);
v_h_acoustic=zeros(1,nb_media.acou);
i_q_acoustic=zeros(1,nb_media.acou);
j_q_acoustic=zeros(1,nb_media.acou);
v_q_acoustic=zeros(1,nb_media.acou);

i_h_eqf=zeros(1,nb_media.eqf);
j_h_eqf=zeros(1,nb_media.eqf);
v_h_eqf=zeros(1,nb_media.eqf);
i_q_eqf=zeros(1,nb_media.eqf);
j_q_eqf=zeros(1,nb_media.eqf);
v_q_eqf=zeros(1,nb_media.eqf);

i_h_limp=zeros(1,nb_media.limp);
j_h_limp=zeros(1,nb_media.limp);
v_h_limp=zeros(1,nb_media.limp);
i_q_limp=zeros(1,nb_media.limp);
j_q_limp=zeros(1,nb_media.limp);
v_q_limp=zeros(1,nb_media.limp);

H_acou=sparse(3*nb_nodes,3*nb_nodes);
Q_acou=sparse(3*nb_nodes,3*nb_nodes);

H_PML=sparse(3*nb_nodes,3*nb_nodes);
Q_PML=sparse(3*nb_nodes,3*nb_nodes);


