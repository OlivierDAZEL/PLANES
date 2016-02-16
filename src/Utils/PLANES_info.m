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
fprintf(fidinfo,'................................................................................\n');
fprintf(fidinfo,'......                           PLANES                                   ......\n');
fprintf(fidinfo,'......                   Porous LAum NumErical Simulator                  ......\n');
fprintf(fidinfo,'................................................................................\n');
fprintf(fidinfo,'......                       Contact: Olivier DAZEL                       ......\n');
fprintf(fidinfo,'......                         LAUM UMR CNRS 6613                         ......\n');
fprintf(fidinfo,'......                     olivier.dazel@univ-lemans.fr                   ......\n');
fprintf(fidinfo,'................................................................................\n');
fprintf(fidinfo,'................................................................................\n');
fprintf(fidinfo,'                                                                                \n');
fprintf(fidinfo,'                                M :  M          8                               \n');
fprintf(fidinfo,'MN.                             $    D Z          M                             \n');
fprintf(fidinfo,'  M.~    M                      NN   M       8  M                               \n');
fprintf(fidinfo,'    7?I      M             M~   M  NM$ MM   M    M                              \n');
fprintf(fidinfo,'       M.+      ,M      M~~                M   7                                \n');
fprintf(fidinfo,'          M  =       M$      :           ZM     M                               \n');
fprintf(fidinfo,'            OI  D               D$        M  D                                  \n');
fprintf(fidinfo,'               $$  8            D  ?MI   M D$    M,                             \n');
fprintf(fidinfo,'                  M,M    M      MO            M ::   M ,                        \n');
fprintf(fidinfo,'                     MM+                    =D     ,M?N                         \n');
fprintf(fidinfo,'                 ~O.N  MM~     M.  Z      N                                     \n');
fprintf(fidinfo,'                Z M   ,    ,M   :     M       ~M8                               \n');
fprintf(fidinfo,'               O   M,           =MM       M          :MD                        \n');
fprintf(fidinfo,'                  +               ZM    NMM     +N          IMD,                \n');
fprintf(fidinfo,'               M  M              MM               MMM7  .:MI        MM,         \n');
fprintf(fidinfo,'                                  8M,M                       $MNM, ,N8. .  MM   \n');
fprintf(fidinfo,'....................................Z...........................................\n');
fprintf(fidinfo,'................................................................................\n');
[~,name_computer]=system('hostname');
fprintf(fidinfo,'...... Generated %s on %s', datestr(now,'dd-mm-yyyy at HH-MM-SS'),name_computer);
fprintf(fidinfo,'...... Name of project = %s\n',project.name);
fprintf(fidinfo,'...... Subproject # %d\n',project.num);
fprintf(fidinfo,'...... #dof FEM= %d\n',nb.dof_FEM);

if nb.dof_FEM~=0
    if (nb.media.elas~=0)
        for i_mat=1:nb.media.elas
            eval(['temp=diag(M_elas_' num2str(i_mat) ');']);
            temp=length(find(temp));
            fprintf(fidinfo,'...... \t\t #dof Elastic solid #%d= %d\n',i_mat,temp);
        end
    end
    if (nb.media.eqf~=0)
        for i_mat=1:nb.media.eqf
            eval(['temp=diag(Q_eqf_' num2str(i_mat) ');']);
            temp=length(find(temp));
            fprintf(fidinfo,'...... \t\t #dof Equivalent fluid #%d= %d\n',i_mat,temp);
        end
    end
    if (nb.media.limp~=0)
        for i_mat=1:nb.media.limp
            eval(['temp=diag(Q_limp_' num2str(i_mat) ');']);
            temp=length(find(temp));
            fprintf(fidinfo,'...... \t\t #dof limp #%d= %d\n',i_mat,temp);
        end
    end
    if (nb.media.pem98~=0)
        for i_mat=1:nb.media.pem98
            eval(['temp=diag(Q_pem98_' num2str(i_mat) '+M_pem98_' num2str(i_mat) ');']);
            temp=length(find(temp));
            fprintf(fidinfo,'...... \t\t #dof PEM98 #%d= %d\n',i_mat,temp);
        end
    end
    if (nb.media.pem01~=0)
        for i_mat=1:nb.media.pem01
            eval(['temp=diag(Q_pem01_' num2str(i_mat) '+M_pem01_' num2str(i_mat) ');']);
            temp=length(find(temp));
            fprintf(fidinfo,'...... \t\t #dof PEM01 #%d= %d\n',i_mat,temp);
        end
    end
    temp=length(find(diag(H_acou)));
    fprintf(fidinfo,'...... \t\t #dof acou= %d\n',temp);
    temp=length(find(diag(H_PML)));
    fprintf(fidinfo,'...... \t\t #dof PML = %d\n',temp);
    fprintf(fidinfo,'...... #dof DGM= %d\n',nb.dof_DGM);
    fprintf(fidinfo,'...... #dof R= %d\n',nb.R);
    fprintf(fidinfo,'...... #dof T= %d\n',nb.T);
    fprintf(fidinfo,'...... #nnz= %d\n',nnz(A));
    fprintf(fidinfo,'...... Computation time = %d s\n',time_PLANES);
    fprintf(fidinfo,'...... Computation time / freq = %d s\n',time_PLANES/length(frequency.vec));
    fprintf(fidinfo,'................................................................................\n');
    fprintf(fidinfo,'................................................................................\n');
    
    fclose(fidinfo);
end