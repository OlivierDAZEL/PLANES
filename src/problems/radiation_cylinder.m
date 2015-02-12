% radiation_cylinder.m
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

lambda=air.c/max(vec_freq);


d_FEM=lambda/12;


r1=0.1;
nr1=ceil(2*pi*r1/d_FEM)
nr1=12;
r2=2.5*r1
nr2=floor(nr1*r2/r1)


dair=3*r1;
dpml=1*dair;
nair=ceil(dair*nr1/(2*pi*r1))
npml=ceil(dpml*nr1/(2*pi*r1))






fid=fopen(name_file_input_FreeFem,'w');
fprintf(fid,'%s\n',name_file_msh);
fprintf(fid,'%12.8f\n',dair);
fprintf(fid,'%12.8f\n',dpml);
fprintf(fid,'%d\n',nair);
fprintf(fid,'%d\n',npml);
fprintf(fid,'%12.8f\n',r1);
fprintf(fid,'%d\n',nr1);
fprintf(fid,'%12.8f\n',r2);
fprintf(fid,'%d\n',nr2);
fclose(fid);

