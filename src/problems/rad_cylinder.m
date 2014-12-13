% rad_cylinder.m
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


lx=0.2;
ly=0.2;
ax=0.2;
ay=0.2;
bx=0.2;
by=0.2;
rayon=0.05;
nlx=2;
nly=2;
nax=2;
nay=2;
nbx=2;
nby=2;
nsrc=16;



fid=fopen(nom_fichier_input_FreeFem,'w');
fprintf(fid,'%s\n',nom_fichier_msh);
fprintf(fid,'%12.8f\n',lx);
fprintf(fid,'%12.8f\n',ly);
fprintf(fid,'%12.8f\n',ax);
fprintf(fid,'%12.8f\n',ay);
fprintf(fid,'%12.8f\n',bx);
fprintf(fid,'%12.8f\n',by);
fprintf(fid,'%12.8f\n',rayon);
fprintf(fid,'%d\n',nlx);
fprintf(fid,'%d\n',nly);
fprintf(fid,'%d\n',nax);
fprintf(fid,'%d\n',nay);
fprintf(fid,'%d\n',nbx);
fprintf(fid,'%d\n',nby);
fprintf(fid,'%d\n',nsrc);
fclose(fid);

