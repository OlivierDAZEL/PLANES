% clear_temporary_FEM_matrices.m
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


clear i_k0_pem01;
clear j_k0_pem01;
clear v_k0_pem01;
clear i_k1_pem01;
clear j_k1_pem01;
clear v_k1_pem01;
clear i_m_pem01;
clear j_m_pem01;
clear v_m_pem01;
clear i_h_pem01;
clear j_h_pem01;
clear v_h_pem01;
clear i_q_pem01;
clear j_q_pem01;
clear v_q_pem01;
clear i_c_pem01;
clear j_c_pem01;
clear v_c_pem01;
clear i_cp_pem01;
clear j_cp_pem01;
clear v_cp_pem01;



clear i_k0_pem98;
clear j_k0_pem98;
clear v_k0_pem98;
clear j_k1_pem98;
clear v_k1_pem98;
clear i_m_pem98;
clear j_m_pem98;
clear v_m_pem98;
clear i_h_pem98;
clear j_h_pem98;
clear v_h_pem98;
clear i_q_pem98;
clear j_q_pem98;
clear v_q_pem98;
clear i_c_pem98;
clear j_c_pem98;
clear v_c_pem98;

clear i_k0_elastic;
clear j_k0_elastic;
clear v_k0_elastic;
clear i_k1_elastic;
clear j_k1_elastic;
clear i_m_elastic;
clear j_m_elastic;
clear v_m_elastic;


clear i_h_acoustic;
clear j_h_acoustic;
clear v_h_acoustic;
clear i_q_acoustic;
clear j_q_acoustic;
clear v_q_acoustic;

clear i_h_eqf;
clear j_h_eqf;
clear v_h_eqf;
clear i_q_eqf;
clear j_q_eqf;
clear v_q_eqf;

clear i_h_limp;
clear j_h_limp;
clear v_h_limp;
clear i_q_limp;
clear j_q_limp;
clear v_q_limp;
