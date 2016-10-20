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


function PLANES(projectname, expnb, data_model, frequency, name)
%
% Runs the PLANES solver for the experiment nb expnb of the
% given project (named projectname, see doc for more info).
%
% projectname -- (string) name of the project
% expnb -- (int) identifier of the experiment (sub-project)
% data_model -- (scalar struct) various information concerning
%                   the model and tools to run
% frequency -- (scalar struct) frequency range of interest
%                   (min, max and number of steps)
% name -- (scalar struct, OPTIONAL ) f specified, name.project_directory must
%                   be a string containing the full path to the project dir

project.name=projectname;
project.num=expnb;

if !isfield(data_model, 'verbosity')
    data_model.verbosity = 0;
end

project.logger = @(level, section, msg) logger(data_model.verbosity, level, section, msg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of PLANES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PLANES_init

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Execution of the data File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

project.logger(3, 'PLANES', 'Running user-defined data script.')
eval([name.project_full '_data']) %

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PLANES_preprocess

if data_model.profiles.mesh
    display_mesh
end

if nb.dof_FEM>0
    EF_global_build
end



PLANES_resolution

PLANES_info

% if (nb.dof_FEM+nb.dof_DGM)>0
%     if (nb.R~=0)
%         fclose(file_abs_id);
%     end
%     if (nb.T~=0)
%         fclose(file_TL_id);
%     end
% end

if exist('FF.inp','file')
    system('rm FF.inp');
end



end

