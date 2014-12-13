% insert_temporary_matrices_elas.m
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


[i_temp,j_temp,s_temp] = find(vk0);
ll=length(i_temp);
i_k0_elas=[i_k0_elas(1:end);index_e(i_temp)'];
j_k0_elas=[j_k0_elas;index_e(j_temp)'];
v_k0_elas=[v_k0_elas;s_temp];
[i_temp,j_temp,s_temp] = find(vk1);
i_k1_elas=[i_k1_elas;index_e(i_temp)'];
j_k1_elas=[j_k1_elas;index_e(j_temp)'];
v_k1_elas=[v_k1_elas;s_temp];
[i_temp,j_temp,s_temp] = find(vm);
i_m_elas=[i_m_elas;index_e(i_temp)'];
j_m_elas=[j_m_elas;index_e(j_temp)'];
v_m_elas=[v_m_elas;s_temp];
