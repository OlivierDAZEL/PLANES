% PLANES_init.m
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

list_path=['''../DGM'',''../FEM'',''../Materials'',''../Mesh'',''../Physics'',''../Plots'',''../Polynomials'',''../Problems'',''../PW'',''../Solutions'',''../Utils'',''../Validation'''];
eval(['addpath(' list_path ');'])


warning('off','MATLAB:nearlySingularMatrix')


name.project_directory=['../../Projects/' project.name '/'];

if project.num==0
    name.file=[name.project_directory project.name];
    name.project_full=project.name;
else
    name.file=[name.project_directory project.name '_' num2str(project.num) ];
    name.project_full=[project.name '_' num2str(project.num)];
end

name.file_msh=          [name.file '.msh'];
name.file_edp=          [name.file '.edp'];
name.file_input_FreeFem=['FF.inp'];
name.file_abs=          [name.file '.abs'];
name.file_L2=           [name.file '.L2'];
name.file_TL=           [name.file '.TL'];
name.file_PW=           [name.file '.PW'];
name.file_info=         [name.file '.info'];
name.file_FEM=          [name.file '.FEM'];

name.compute_error=['compute_error_' , project.name];
if ((exist(name.compute_error)==2))
    name.file_error=          [name.file '.error'];
end

name.solution=['solution_' , project.name];


if export.profiles==1
    name.directory_profiles= [name.project_directory '/Profiles/'];
    if ~exist(name.directory_profiles,'dir')
        mkdir(name.directory_profiles)
    end
end

set(0,'DefaultLineMarkerSize',15);
set(0,'Defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',15)
set(0,'DefaultAxeslinewidth',2);




Gauss_points=compute_GaussLegendre_points(10);
