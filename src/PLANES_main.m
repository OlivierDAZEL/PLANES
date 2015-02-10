% PLANES_main.m
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

clear all
close all
clc
list_path=['''FEM'',''DGM'',''problems'',''Materials'',''Mesh'',''Physics'',''plots'',''Utils'',''validation'',''PW'',''analytical_solutions'''];
eval(['addpath(' list_path ');'])

name_of_project='air_PEM';
subproject=0;
% Number of frequencies
% If the number is negative then a logscale is chosen
% If the number is equal to 1, then the frequency is equal to freq_min
nb_frequencies=1;
freq_min=100;
freq_max=4000;
% Angle of incidence
theta=0*pi/180;

solve.FEM=1;
solve.DGM=0;
solve.PW=0;

if solve.DGM
    nb_theta=4;
end

%
export_profiles=0;
plot_profiles=1;
export_nrj=0;

%% Initialization of PLANES
init_PLANES
air_properties_JPG
init_vec_frequencies


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call to FreeFEM++ to create the Mesh and importation of the mesh
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eval([name_of_project_full])

system(['/usr/local/bin/FreeFem++ ' name_file_edp])

tic

disp('Importing Mesh')

if solve.DGM
    
    [nb_nodes,nb_elements,nb_edges,nodes,elements,element_label,edges,nb_media,num_media,element_num_mat,nb_internal,internal,...
        nb_MMT,edges_MMT,nb_loads,loads,nb_dirichlets,dirichlets,nb_periodicity,periodicity]=msh2DGM(name_file_msh,1);
    
    analyze_mesh_DGM
    
end

if solve.FEM
    
    [nb_nodes,nb_elements,nb_edges,nodes,elements,element_label,edges,nb_media,num_media,element_num_mat,nb_interfaces,interfaces,...
        nb_MMT,edges_MMT,nb_loads,loads,nb_dirichlets,dirichlets,nb_periodicity,periodicity]=msh2TR6(name_file_msh,0);
    analyze_mesh_FEM
    disp('Building FEM shape matrices')
    
    EF_TR6_global_build
    
end


disp('End of mesh importation')
toc

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction of global shape matrices
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency loop
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if solve.FEM
    disp('FEM Resolution')
    FEM_resolution
end

if solve.DGM
    disp('DGM Resolution')
    DGM_resolution
end

if solve.PW
    disp('PW Resolution')
    PW_resolution
end

% plot_sol_1D_r
% 
% Amp_cylinder=(2*air.rho*omega^2)/(k_air*(besselh(-1,2,k_air*r1)-besselh(1,2,k_air*r1)));
% r_plot=linspace(r1,r2,2000);
% p_plot=Amp_cylinder*besselh(0,2,k_air*r_plot);
% plot(r_plot,abs(p_plot),'k')



 solution_AIR_PEM





% A=1/(j*omega*sin(k_air*ly));
% x_analytique=linspace(-ly,0,200);
% sol_analytique=-air.K*k_air*A*cos(k_air*x_analytique);
%figure(10002)
%plot(x_analytique+ly,abs(sol_analytique))

% figure
% hold on
% plot(vec_freq,abs_EF,'m.')
%plot(vec_freq,abs_PW)

% figure
% semilogx(vec_freq,TL_EF,'.')
% hold on
% semilogx(vec_freq,TL_PW)

%eval(['rmpath(' list_path ');'])