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


%function PLANES_Multilayer_3D(projectname, expnb,data_model,multilayer,frequency)
clear all;
close all;
clc;



frequency.nb=-100;
frequency.min=1;
frequency.max=20000;


multilayer(1,1).nb=1;
multilayer(1,1).d=5e-2;
multilayer(1,1).mat=5010;
multilayer(1,1).termination=0;

data_model.theta_1=45*pi/180;
data_model.theta_2=62*pi/180;




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
    project.num=101;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of PLANES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PLANES_init
PLANES_preprocess


PW_3D_resolution

PLANES_info

figure 
semilogx(frequency.vec,abs_PW_3D,'.-')
load('../../Projects/Multilayer_3D/out/Multilayer_3D_101.PW3D')
load('../../Projects/Multilayer_3D/test.mat');
hold on
semilogx(test(1).f,test(1).absorp,'r')

maine=load('../../../../Programmation/Maine/TCLTK/out.dat')
semilogx(maine(:,1),maine(:,4),'k')



%end
