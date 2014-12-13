% Benchmark6.m
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


theta=60*pi/180;
dxair=2e-2;
dyair=5.5e-2;
nx=10;
ny=10;
labelexcitation=12;
labelporous=5004;
r=2.5e-3;
labelelastic=1001;
ninclusion=2*nx;


fid=fopen(name_file_input_FreeFem,'w');
fprintf(fid,'%s\n',name_file_msh);
fprintf(fid,'%12.8f\n',dxair);
fprintf(fid,'%12.8f\n',dyair);
fprintf(fid,'%d\n',nx);
fprintf(fid,'%d\n',ny);
fprintf(fid,'%d\n',labelexcitation);
fprintf(fid,'%d\n',labelporous);
fprintf(fid,'%12.8f\n',r);
fprintf(fid,'%d\n',labelelastic);
fprintf(fid,'%d\n',ninclusion);
fclose(fid);