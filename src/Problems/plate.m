% sandwich_meta.m
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
warning('off','all');

period=1e-2;
thicknessplate=1e-3;
labelplate=1001;

nplate=2;
nx=ceil(nplate*period/thicknessplate);

fid=fopen(name.file_input_FreeFem,'w');
fprintf(fid,'%s\n',name.file_msh);
fprintf(fid,'%12.8f\n',period);
fprintf(fid,'%12.8f\n',thicknessplate);
fprintf(fid,'%d\n',labelplate);
fprintf(fid,'%d\n',nx);
fprintf(fid,'%d\n',nplate);
fclose(fid);

% Call to FreeFem++ to create a msh File
system(['/usr/local/bin/FreeFem++ ' name.file_edp]);

%Importation of the msh File

[nb,nodes,elem,edge_msh]=msh_import(name.file_msh);

% All the elements are TR6
elem.model=1*ones(nb.elements,1);

theta_inc=60*pi/180;

% Mat_elas_1
% incident(1).typ=1;
% incident(1).lambda=lambda_solide;
% incident(1).mu=mu_solide;
% incident(1).rho=rho_solide;
% incident(1).thickness=thicknessplate;
% 
% transmitted(1).typ=1;
% transmitted(1).lambda=lambda_solide;
% transmitted(1).mu=mu_solide;
% transmitted(1).rho=rho_solide;
% transmitted(1).thickness=thicknessplate;

nb_layers=1;
multilayer(1).d=thicknessplate;
multilayer(1).mat=labelplate;
% Termination condition // 0 for rigid backing 1 for radiation
termination=1;

% Number of waves (two ways) in each layer
% For the resolution, the incident waves is included in the system and put to RHS at the
% end of the procedure

% Addition of a new layer for the incident medium
l0.d=0;
l0.mat=0;
multilayer=[l0 multilayer];
nb_layers=nb_layers+1;

compute_number_PW_TMM


