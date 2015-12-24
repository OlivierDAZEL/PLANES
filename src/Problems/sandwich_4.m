% sandwich_4.m
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


lx=10e-2;
thicknessplate=1e-3;
thicknessporous=2e-2;
labelplate=1001;
labelporous=5003;
labelbottom=2100;
labeltop=2101;

hbottom=3e-2;
htop=3e-2;


nplate=1;
nporous=1%ceil(thicknessporous*nplate/thicknessplate);
nx=ceil(lx*nporous/thicknessporous);
nbottom=ceil(hbottom*nx/lx);
ntop=ceil(htop*nx/lx);


fid=fopen(name.file_input_FreeFem,'w');
fprintf(fid,'%s\n',name.file_msh);
fprintf(fid,'%12.8f\n',lx);
fprintf(fid,'%12.8f\n',thicknessplate);
fprintf(fid,'%12.8f\n',thicknessporous);
fprintf(fid,'%d\n',labelplate);
fprintf(fid,'%d\n',labelporous);
fprintf(fid,'%d\n',labelbottom);
fprintf(fid,'%d\n',labeltop);
fprintf(fid,'%12.8f\n',hbottom);
fprintf(fid,'%12.8f\n',htop);
fprintf(fid,'%d\n',nplate);
fprintf(fid,'%d\n',nporous);
fprintf(fid,'%d\n',nx);
fprintf(fid,'%d\n',nbottom);
fprintf(fid,'%d\n',ntop);
fclose(fid);

% Call to FreeFem++ to create a msh File
system(['/usr/local/bin/FreeFem++ ' name.file_edp]);

%Importation of the msh File

[nb,nodes,elem,edge_msh]=msh_import(name.file_msh);

% All the elements are TR6
elem.model=1*ones(nb.elements,1);

theta_inc=40*pi/180;


nb_layers=3;
multilayer(1).d=thicknessplate;
multilayer(1).mat=labelplate;
multilayer(2).d=thicknessporous;
multilayer(2).mat=labelporous;
multilayer(3).d=thicknessplate;
multilayer(3).mat=labelplate;
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