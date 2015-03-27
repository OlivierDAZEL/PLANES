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
%close all
clc
list_path=['''FEM'',''DGM'',''problems'',''Materials'',''Mesh'',''Physics'',''plots'',''Utils'',''validation'',''PW'',''analytical_solutions'''];
eval(['addpath(' list_path ');'])

abs_JPG=load('Abs_d55_angpi90_type3.dat');


%name_of_project='IJNME_DGM';
%name_of_project='PEM';
%name_of_project='air_PEM';
%name_of_project='Kundt';
name_of_project='TW';


subproject=5;
% Number of frequencies
% If the number is negative then a logscale is chosen
% If the number is equal to 1, then the frequency is equal to freq_min
nb_frequencies=200;
freq_min=10;
freq_max=5000;
% Angle of incidence
theta_inc=0*pi/180;

solve.FEM=1;
solve.DGM=0;
solve.PW=0;

if solve.DGM
    nb_theta=32;
end

%
export_profiles=1;
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
        nb_MMT,edges_MMT,nb_loads,loads,nb_dirichlets,dirichlets,nb_periodicity,periodicity,index_label,index_element]=msh2DGM(name_file_msh,0);
    
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

%trace_DGM_y

% Analytical solution (if exists)

name_solution=['solution_' , name_of_project];
if ((exist(name_solution)==2)&(plot_profiles))
    eval('eval(name_solution)')
end




open('essai.fig')
hold on
plot(vec_freq,abs_EF,'k.-')

% 
% open('abs_coque_corr.fig')
% hold on
% plot(vec_freq,abs_EF,'.')

% X(end)
% Reflexion_coefficient
%
% X(end)-Reflexion_coefficient

%
% disp('Modules')
% disp('DGM')
% abs(reflex_DGM)
% disp('PW')
% abs(reflex_PW)
%
% disp('Angles')
% disp('DGM')
% angle(reflex_DGM)
% disp('PW')
% angle(reflex_PW)
%
% %
% error_DGM=abs((X(end)-reflex_PW)/X(end))
%
% sol_maine=load('../../../Programmation/Maine/TCLTK/out.dat');

% 
% if solve.FEM
% sol=zeros(3*nb_nodes);
% sol(dof_back)=X(1:end-1);
% save('-v7.3','../../../../Desktop/additional_example/DGM_additional_example_FEM','nodes','elements','sol','reflex_FEM','reflex_PW','nb_dof_FEM','omega')
% end
% 
% if solve.DGM
% save('-v7.3','../../../../Desktop/additional_example/DGM_additional_example_DGM','nodes','elements','X','reflex_DGM','reflex_PW','nb_dof_DGM',...
%     'delta_1','delta_2','delta_3','mu_1','mu_2','mu_3','N','A_hat','K_eq_til','omega','dof_start_element','ondes_element','nb_theta','vec_theta')
% end