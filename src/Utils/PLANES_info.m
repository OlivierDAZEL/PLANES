% PLANES_info.m
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


fidinfo=fopen(name.file_info,'w');
fprintf(fidinfo,'................................................................................');
fprintf(fidinfo,'......                           PLANES                                   ......');
fprintf(fidinfo,'......                   Porous LAum NumErical Simulator                  ......');
fprintf(fidinfo,'................................................................................');
fprintf(fidinfo,'......                       Contact: Olivier DAZEL                       ......');
fprintf(fidinfo,'......                         LAUM UMR CNRS 6613                         ......');
fprintf(fidinfo,'......                     olivier.dazel@univ-lemans.fr                   ......');
fprintf(fidinfo,'................................................................................');
fprintf(fidinfo,'................................................................................');
fprintf(fidinfo,'                                                                                ');
fprintf(fidinfo,'                                M :  M          8                               ');
fprintf(fidinfo,'MN.                             $    D Z          M                             ');
fprintf(fidinfo,'  M.~    M                      NN   M       8  M                               ');
fprintf(fidinfo,'    7?I      M             M~   M  NM$ MM   M    M                              ');
fprintf(fidinfo,'       M.+      ,M      M~~                M   7                                ');
fprintf(fidinfo,'          M  =       M$      :           ZM     M                               ');
fprintf(fidinfo,'            OI  D               D$        M  D                                  ');
fprintf(fidinfo,'               $$  8            D  ?MI   M D$    M,                             ');
fprintf(fidinfo,'                  M,M    M      MO            M ::   M ,                        ');
fprintf(fidinfo,'                     MM+                    =D     ,M?N                         ');
fprintf(fidinfo,'                 ~O.N  MM~     M.  Z      N                                     ');
fprintf(fidinfo,'                Z M   ,    ,M   :     M       ~M8                               ');
fprintf(fidinfo,'               O   M,           =MM       M          :MD                        ');
fprintf(fidinfo,'                  +               ZM    NMM     +N          IMD,                ');
fprintf(fidinfo,'               M  M              MM               MMM7  .:MI        MM,         ');
fprintf(fidinfo,'                                  8M,M                       $MNM, ,N8. .  MM   ');
fprintf(fidinfo,'....................................Z...........................................');
fprintf(fidinfo,'................................................................................');
[~,name_computer]=system('hostname');
fprintf(fidinfo,'...... Generated %s on %s', datestr(now,'dd-mm-yyyy at HH-MM-SS'),name_computer);
fprintf(fidinfo,'...... Name of project = %s\n',project.name);
fprintf(fidinfo,'...... Subproject # %d\n',project.num);
fprintf(fidinfo,'...... #dof FEM= %d\n',nb.dof_FEM);
fprintf(fidinfo,'...... #dof DGM= %d\n',nb.dof_DGM);
fprintf(fidinfo,'...... #dof R= %d\n',nb.R);
fprintf(fidinfo,'...... #dof T= %d\n',nb.T);
fprintf(fidinfo,'...... Computation time = %d s\n',time_PLANES);
fprintf(fidinfo,'....................................Z...........................................');
fprintf(fidinfo,'................................................................................');

fclose(fidinfo);