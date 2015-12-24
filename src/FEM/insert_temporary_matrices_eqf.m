% insert_temporary_matrices_eqf.m
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


[i_temp,j_temp,s_temp] = find(vh);
ll=length(i_temp);
i_h_eqf=[i_h_eqf;index_p(i_temp)'];
j_h_eqf=[j_h_eqf;index_p(j_temp)'];
v_h_eqf(1:end+ll,elem.num_mat(ie))=[v_h_eqf(:,elem.num_mat(ie));s_temp];
[i_temp,j_temp,s_temp] = find(vq);
ll=length(i_temp);
i_q_eqf=[i_q_eqf;index_p(i_temp)'];
j_q_eqf=[j_q_eqf;index_p(j_temp)'];
v_q_eqf(1:end+ll,elem.num_mat(ie))=[v_q_eqf(:,elem.num_mat(ie));s_temp];

