% periodicity_condition_application.m
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


%disp('Applying periodicity condtions')

edge_left= find(edges.periodicity(:,4)==98);
edge_right=find(edges.periodicity(:,4)==99);

node_left=unique([edges.periodicity(edge_left,1);edges.periodicity(edge_left,2);edges.periodicity(edge_left,6)]);
[temp,i_left]=sort(nodes(node_left,2));
node_left=node_left(i_left);

node_right=unique([edges.periodicity(edge_right,1);edges.periodicity(edge_right,2);edges.periodicity(edge_right,6)]);
[temp,i_right]=sort(nodes(node_right,2));
node_right=node_right(i_right);

dof_left= dof_A(uxyp_TR(node_left));
dof_right=dof_A(uxyp_TR(node_right));

dof_left= dof_left (find(dof_left));
dof_right=dof_right(find(dof_right));

dof_internal=1:nb.dof_total;
dof_internal([dof_left dof_right])=[];


AA=A(dof_left,dof_left);
BB=A(dof_left,dof_internal);
CC=A(dof_left,dof_right);
DD=A(dof_internal,dof_left);
EE=A(dof_internal,dof_internal);
FF=A(dof_internal,dof_right);
GG=A(dof_right,dof_left);
HH=A(dof_right,dof_internal);
II=A(dof_right,dof_right);

Mat_M=[CC 0*HH;0*FF 0*EE]; % delta^2
Mat_C=[AA+II BB;FF 0*EE];  %delta
Mat_K=[GG HH;DD EE];

AAA=[Mat_M 0*Mat_M;0*Mat_M -Mat_K];
BBB=[Mat_C   Mat_K; Mat_K 0*Mat_K];


delta=eigs(inv(BBB)*AAA,10,'SR')

k=log(delta)/(-1j*data_model.L)



fsddfdsfdsdsf








