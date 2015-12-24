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

% project.name='Kundt';
% project.num=1;
% project.name='sandwich_meta';
% project.name='plate';
% project.name='sandwich';
project.name='Article_FEM_TMM';
project.num=1;

%project.name='sandwich_meta';
%project.name='sandwich_siniat';
%project.num=3;

% project.name='Incompatible_mesh';
%project.name='CFM';
%project.name='CFM_refined';
%project.num=1;

% Number of frequencies
% If the number is negative then a logscale is chosen
% If the number is equal to 1, then the frequency is equal to freq_min

frequency.nb=1;
frequency.min=10;
frequency.max=1000;

profiles.mesh=0;
profiles.solution=0;
profiles.x=0;
profiles.y=0;
profiles.map=0;
profiles.custom=0;
profiles.on=profiles.x+profiles.y+profiles.map+profiles.custom;
export.nrj=1;
export.profiles=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of PLANES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PLANES_init
air_properties_maine
init_vec_frequencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creation and importation of the mesh
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eval([name.project_full]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PLANES_preprocess

if nb.dof_FEM>0
    EF_global_build
end

PLANES_resolution

PLANES_info

%Analytical solution (if exists)
if ((exist(name.solution)==2)&&(profiles.on~=0))
    eval('eval(name.solution)')
end

%Maine=load('../../../../Programmation/Maine/TCLTK/out.dat');

% figure(1)
% %semilogx(Maine(:,1),Maine(:,2))
% %hold on
% semilogx(frequency.vec,TL_PW,'r.-')
% hold on
% if project.num==0
%     semilogx(frequency.vec,TL_EF,'k')
%     semilogx(frequency.vec,abs(TL_EF-TL_PW),'k')
% else
%     semilogx(frequency.vec,TL_EF,'m')
%     semilogx(frequency.vec,abs(TL_EF-TL_PW),'m')
% end


figure(666666)

 if project.num==10
     plot(frequency.vec,Db,'k')
     hold on

 else
     plot(frequency.vec,Db,'m.-')
     hold on
 end
 
%  TL_PW
%  TL_EF
%  abs(TL_PW-TL_EF)
%  
%  Db(1)
% target=0.916939036846682

% figure(2010)
% k_eq=omega*sqrt(rho_eq_til/K_eq_til);
% x_m=thicknessplate1+thicknessplate2+thicknessporous;
% Zs=-j*sqrt(rho_eq_til*K_eq_til)*cot(k_eq*x_m);
% R=(Zs-air.Z)/(Zs+air.Z)
% rflx
% rflx_PW
% rflx_PW-rflx

%A=-1/sin(k_eq*(x_m))/(j*omega);

% x=linspace(0,x_m,500);
% 
% p=1*exp(-1j*k_eq*x)+R*exp(1j*k_eq*x);


%plot(x,abs(p),'k')
