% File PLANES_visualize_profiles.m
%
% Copyright (C) 2014 Olivier DAZEL <olivier.dazel@univ-lemans.fr>
%
% This file is part of PLANES (Porous LAum Numerical Simulator)  
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
% 
% More details can be found on http://perso.univ-lemans.fr/~odazel/
%
% Up-to-date copies of PLANES can be obtained from the Web page:
% https://github.com/OlivierDAZEL/PLANES


% Initialization of PLANES

clear all
close all
clc
list_path=['''FEM'',''generate_mesh'',''Materials'',''Mesh'',''Physics'',''plots'',''Utils'',''validation'''];
eval(['addpath(' list_path ');'])

name_of_project='Benchmark1';
subproject=0;


name_project_directory=['../Projects/' name_of_project '/'];
name_directory_profiles= [name_project_directory '/Profiles/'];
load([name_directory_profiles 'vec_freq.mat'])
load([name_directory_profiles 'elements.mat'])
load([name_directory_profiles 'element_label.mat'])
nb_elements=size(elements,1);
load([name_directory_profiles 'nodes.mat'])
load([name_directory_profiles 'pressures.mat'])

max_pressure=max(max(abs(pressures)));


for i_f=1:1%length(vec_freq)   
    figure
    set(gca,'Fontsize',20)
    hold on
    sol=pressures(:,i_f);
    plot_solution_TR6
    axis equal
    caxis([0 max_pressure])
    set(gca,'xtick',[min(nodes(:,1)) max(nodes(:,1))]);
    set(gca,'ytick',[min(nodes(:,2)) max(nodes(:,2))]);
    colorbar('Fontsize',25) 
    title(['Pressure at ' num2str(vec_freq(i_f)) ' Hz'],'Fontsize',20)

    print('-djpeg99',[name_directory_profiles num2str(i_f)])
    %close
end





%eval(['rmpath(' list_path ');'])