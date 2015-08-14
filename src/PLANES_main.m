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

%name_of_project='PEM';
%name_of_project='air_PEM';
project.name='Kundt';
%name_of_project='CFM';
%name_of_project='sandwich_meta';
project.num=4;

%subproject=1;
% Number of frequencies
% If the number is negative then a logscale is chosen
% If the number is equal to 1, then the frequency is equal to freq_min
frequency.nb=1;
frequency.min=2000;
frequency.max=4000;
theta_DGM.nb=2;
vec_theta=linspace(0,2*pi,theta_DGM.nb+1);
vec_theta(end)=[];
vec_theta=vec_theta+pi/2;

% Angle of incidence
%theta_inc=0*pi/180;



% if (solve.DGM|solve.FEMDGM)
%     nb_thetaDGM=8;
%     vec_theta=linspace(0,2*pi,nb_thetaDGM+1);
%     vec_theta(end)=[];
% end
%


profiles.mesh=1;
profiles.x=0;
profiles.y=1;
profiles.map=0;
profiles.custom=0;
profiles.on=profiles.x+profiles.y+profiles.map+profiles.custom;
export.nrj=1;
export.profiles=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of PLANES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PLANES_init
air_properties_JPG
init_vec_frequencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call to FreeFEM++ to create the Mesh and importation of the mesh
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eval([name.project_full]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



PLANES_preprocess




if profiles.mesh
    display_mesh
end




if nb.dof_FEM>0
    EF_global_build
end


PLANES_resolution





if solve.FEMDGM
    
    [nb,nodes,elements,element_label,element_model,edges,num_media,element_num_mat,internal,...
        edges_MMT,loads,dirichlets,periodicity,index_label,index_element,lx_H12,ly_H12]=createmshH12DGM(lx,ly/2,ly/2,nx,ny/2,ny/2,1);
    %element_model=ones(size(element_model));
    
    analyze_mesh_H12
    analyze_mesh_DGM
    
    EF_H12_global_build
    

    FEMDGM_resolution
 
end


if solve.PW
    disp('PW Resolution')
     PW_resolution
end


% Analytical solution (if exists)

name_solution=['solution_' , project.name];
if project.num~=0
    name_solution=[name_solution,'_',num2str(project.num)];
end
if ((exist(name_solution)==2)&(profiles.on))
    eval('eval(name_solution)')
end

