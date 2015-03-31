% plot_sol_TR6_custom_55.m
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



x=[];
y=[];


center_inclusion=[period/2 thicknessporous/2];


figure(10055)
subplot(221)
hold on
subplot(222)
hold on
subplot(223)
hold on
subplot(224)
hold on


for ie=1:nb.interfaces
    node=interfaces(ie,1);
    theta_1=angle((nodes(node,1)-center_inclusion(1))+1j*(nodes(node,2)-center_inclusion(2)));
    
    
    
    ux_1=sol(ux(node));
    uy_1=sol(uy(node));
    u_1=sqrt(abs(sol(ux(node)))^2+abs(sol(uy(node)))^2);
    ur_1=[ux_1 uy_1]*[cos(theta_1);sin(theta_1)];
    utheta_1=[ux_1 uy_1]*[-sin(theta_1);cos(theta_1)];
    node=interfaces(ie,2);
    theta_2=angle((nodes(node,1)-center_inclusion(1))+1j*(nodes(node,2)-center_inclusion(2)));
    
    
    
    ux_2=sol(ux(node));
    uy_2=sol(uy(node));
    u_2=sqrt(abs(sol(ux(node)))^2+abs(sol(uy(node)))^2);
    ur_2=[ux_2 uy_2]*[cos(theta_2);sin(theta_2)];
    utheta_2=[ux_2 uy_2]*[-sin(theta_2);cos(theta_2)];
    
    node=interfaces(ie,6);
    
    
    theta_3=angle((nodes(node,1)-center_inclusion(1))+1j*(nodes(node,2)-center_inclusion(2)));
    ux_3=sol(ux(node));
    uy_3=sol(uy(node));
    u_3=sqrt(abs(sol(ux(node)))^2+abs(sol(uy(node)))^2);
    ur_3=[ux_3 uy_3]*[cos(theta_3);sin(theta_3)];
    utheta_3=[ux_3 uy_3]*[-sin(theta_3);cos(theta_3)];
    
    figure(10055)
    subplot(221)
    polar([theta_1 theta_3 theta_2],abs([u_1 u_3 u_2]),'k.')
    subplot(222)
    polar([theta_1 theta_3 theta_2],abs([ux_1 ux_3 ux_2]),'r.')
    polar([theta_1 theta_3 theta_2],abs([uy_1 uy_3 uy_2]),'m.')
    
    
    subplot(223)
    polar([theta_1 theta_3 theta_2],abs([ur_1 ur_3 ur_2]),'r.')
    polar([theta_1 theta_3 theta_2],abs([utheta_1 utheta_3 utheta_2]),'m.')
    subplot(224)
    plot(xd,yd,'m')
    plot(vec_freq(i_f),abs_EF(i_f),'k.')
    
end

 subplot(221)
 axis equal
subplot(222)
 axis equal
subplot(223)
 axis equal
 


if export_profiles==1
    print('-djpeg',[name_directory_profiles, num2str(i_f)]);
end

close(10055)
