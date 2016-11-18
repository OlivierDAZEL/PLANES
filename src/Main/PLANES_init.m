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

list_path=['''../DGM'',''../FEM'',''../Mesh'',''../Physics'',''../Plots'',''../Polynomials'',''../PW'',''../Solutions'',''../Utils'''];
eval(['addpath(' list_path ');'])


warning('off','all');
%warning('off','MATLAB:nearlySingularMatrix')

if exist('name')~=1
    name.project_directory=['../../Projects/' project.name '/'];
end
addpath([name.project_directory '/m/']);
addpath([name.project_directory '/m/Materials']);
if name.project_directory(end)~='/'
    name.project_directory = [name.project_directory '/'];
end

if (project.num)==0
    name.project_full=project.name;
else
    name.project_full=[project.name '_' num2str(project.num)];
end


name.file_msh=          [name.project_directory 'FF/'   name.project_full '.msh'];
name.file_edp=          [name.project_directory 'FF/'   name.project_full '.edp'];
name.file_FF=           [name.project_directory 'FF/FF.inp'];
name.file_abs=          [name.project_directory 'out/'  name.project_full '.abs'];
name.file_L2=           [name.project_directory 'out/'  name.project_full '.L2'];
name.file_PL=           [name.project_directory 'out/'  name.project_full '.PL'];
name.file_TL=           [name.project_directory 'out/'  name.project_full '.TL'];
name.file_PW=           [name.project_directory 'out/'  name.project_full '.PW'];
name.file_info=         [name.project_directory 'out/'  name.project_full '.info'];
name.file_FEM=          [name.project_directory 'out/'  name.project_full '.FEM'];
name.dir_out=           [name.project_directory 'out/'];

set(0,'DefaultLineMarkerSize',15);
set(0,'Defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',15)
set(0,'DefaultAxeslinewidth',2);

eval(['air_properties_' data_model.air_model ])

if ~isfield(data_model, 'ngauss')
	data_model.ngauss = 10;
end
Gauss_points=compute_GaussLegendre_points(data_model.ngauss);
project.logger(1, 'init', ['Gauss quadrature using ' num2str(data_model.ngauss) ' points']);

if data_model.export.reset~=0;
   system(['rm ' name.file_PL ';']); 
end
