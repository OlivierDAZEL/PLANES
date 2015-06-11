% init_PLANES.m
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


warning('off','MATLAB:nearlySingularMatrix')


name_project_directory=['../Projects/' name_of_project '/'];

if subproject==0
    name_file=[name_project_directory name_of_project];
    name_of_project_full=name_of_project;
else
    name_file=[name_project_directory name_of_project '_' num2str(subproject) ];
    name_of_project_full=[name_of_project '_' num2str(subproject)];
end

name_file_msh=          [name_file '.msh'];
name_file_edp=          [name_file '.edp'];
name_file_input_FreeFem=['FF.inp'];
name_file_abs=          [name_file '.abs'];
name_file_TL=           [name_file '.TL'];
name_file_PW=           [name_file '.PW'];
name_file_info=         [name_file '.info'];

name_file_FEM=          [name_file '.FEM'];

if solve.DGM
    name_file_DGM=          [name_file '_Nw=' num2str(nb_thetaDGM) '.DGM'];
end

if export_profiles==1
    name_directory_profiles= [name_project_directory '/Profiles/'];
    if ~exist(name_directory_profiles,'dir')
        mkdir(name_directory_profiles)
    end
end

set(0,'DefaultLineMarkerSize',15);
set(0,'Defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',15)
set(0,'DefaultAxeslinewidth',2);