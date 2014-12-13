% PLANES_Benchmark.m
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


clear all
close all
clc
cd('Utils')
    init_path
cd ..

export_profiles=1;
nb_frequencies=1;


name_of_project='Benchmark1';
subproject=0;
plot_abs=1;
plot_TL=0;
init_PLANES

if plot_abs==1
    PLANES_result=load([name_project_directory name_of_project '.abs']);
    reference=load([name_project_directory name_of_project '.ref']);
    figure
    hold on
    plot(PLANES_result(:,1),PLANES_result(:,2),'r.')
    plot(reference(:,1),reference(:,4))
end
    
    