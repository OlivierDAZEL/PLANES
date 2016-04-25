% plot_PLANES_map.m
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


if sum(ismember(floor(elem.label/1000),[1 4 5]))~=0
    figure(10001)
    hold on
    title('Solid displacement')
    xlabel('x')
    ylabel('y')
%     figure(11002)
%     hold on
%     title('Solid displacement')
%     xlabel('x')
%     ylabel('y')
end


if sum(ismember(floor(elem.label/1000),[0 2 3 4 5]))~=0
    figure(10002)
    hold on
    title('Pressure')
    xlabel('x')
    ylabel('y')
%     figure(11010)
%     hold on
%     title('Pressure')
%     xlabel('y')
%     ylabel('angle(P)')
end



for ie=1:nb.elements
    if ismember(elem.model(ie),[1 4])
            plot_sol_TR6_on_element
   elseif ismember(elem.model(ie),[10 11])
            plot_sol_DGM_on_element
    elseif ismember(elem.model(ie),[2])
            plot_sol_H12_on_element
    end
end



if data_model.export.profiles==1
    figure(10002)    
    shading interp
    colormap jet
    colorbar
    print('-djpeg',[name.directory_profiles, num2str(i_f) ,'without_mesh']);
    display_mesh_light
    print('-djpeg',[name.directory_profiles, num2str(i_f) ,'with_mesh']);
    close(10002) 
    figure(10002)    
    display_mesh_light
    print('-djpeg',[name.directory_profiles, num2str(i_f) ,'only_mesh']);

end

if sum(ismember(floor(elem.label/1000),[1 4 5]))~=0
    figure(10001)
    colormap jet
    colorbar
end


if sum(ismember(floor(elem.label/1000),[0 2 3 4 5]))~=0
    figure(10002)
    colormap jet
    colorbar
%     figure(11002)
%     colormap jet
%     colorbar
end

