% PLANES.m
%
% Copyright (C) 2016 < Olivier DAZEL <olivier.dazel@univ-lemans.fr> >
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
%                 https://github.com/OlivierDAZEL/PLANES
% or find more details on Olivier's webpage
%                  http://perso.univ-lemans.fr/~odazel/
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
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%


function PLANES(projectname, expnb,data_model,frequency,profiles,export)
%clear all;
%close all;
%clc;

if ~exist('projectname')==1
    project.name='David_ZOD';
    project.name='Article_ZOD';
    project.name='Multilayer_3D';
    nargin=0;
else
    project.name=projectname;
end
if exist('expnb')==1
    project.num=expnb;
else
    project.num=21;
    project.num=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of PLANES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PLANES_init

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Execution of the data File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eval([name.project_full '_data']) %



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


PLANES_preprocess


if profiles.mesh
    display_mesh
end

if nb.dof_FEM>0
    EF_global_build
end

PLANES_resolution

PLANES_info

if (nb.dof_FEM+nb.dof_DGM)>0
    if (nb.R~=0)
        fclose(file_abs_id);
    end
    if (nb.T~=0)
        fclose(file_TL_id);
    end
end

if exist('FF.inp','file')
    system('rm FF.inp')
end


%end

maine=load('../../../../Programmation/Maine/TCLTK/out.dat');
figure
semilogx(maine(:,1),maine(:,4))
hold on
semilogx(frequency.vec,abs_PW,'r.')
semilogx(frequency.vec,abs_PW_3D,'r')
end

% load('OD1')
% semilogx(f,1-abs(R).^2,'k--')


% figure
% semilogx(frequency.vec,abs(abs_PW_3D-maine(:,4)'),'r')
% hold on
%semilogx(f,abs(1-abs(R).^2-maine(:,4)'),'k')


