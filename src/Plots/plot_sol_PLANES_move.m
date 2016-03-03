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
factor=0.05/max(abs(sol))

if sum(ismember(floor(elem.label/1000),[1 4 5]))~=0
    figure(5000)
    hold on
    title('Solid displacement')
    xlabel('x')
    ylabel('y')
end





for ie=1:nb.elements
    x=nodes(elem.nodes(ie,1:6),1);
    y=nodes(elem.nodes(ie,1:6),2);

    
    if elem.model(ie)==1
        line([x(1) x(3)],[y(1) y(3)],'Color','g','Linewidth',0.5);
        line([x(3) x(5)],[y(3) y(5)],'Color','g','Linewidth',0.5);
        line([x(5) x(1)],[y(5) y(1)],'Color','g','Linewidth',0.5);
    end
   
    x_move=factor*(imag(sol(ux_TR(elem.nodes(ie,1:6)))));
    y_move=factor*(imag(sol(uy_TR(elem.nodes(ie,1:6)))));
    
    
    if elem.model(ie)==1
        line([x(1)+x_move(1) x(2)+x_move(2)],[y(1)+y_move(1) y(2)+y_move(2)],'Color','r','Linewidth',2);
        line([x(2)+x_move(2) x(3)+x_move(3)],[y(2)+y_move(2) y(3)+y_move(3)],'Color','r','Linewidth',2);
        line([x(3)+x_move(3) x(4)+x_move(4)],[y(3)+y_move(3) y(4)+y_move(4)],'Color','r','Linewidth',2);
        line([x(4)+x_move(4) x(5)+x_move(5)],[y(4)+y_move(4) y(5)+y_move(5)],'Color','r','Linewidth',2);
        line([x(5)+x_move(5) x(6)+x_move(6)],[y(5)+y_move(5) y(6)+y_move(6)],'Color','r','Linewidth',2);
        line([x(6)+x_move(6) x(1)+x_move(1)],[y(6)+y_move(6) y(1)+y_move(1)],'Color','r','Linewidth',2);
      
    end
   
    x_move=factor*(real(sol(ux_TR(elem.nodes(ie,1:6)))));
    y_move=factor*(real(sol(uy_TR(elem.nodes(ie,1:6)))));
    
    
    if elem.model(ie)==1
        line([x(1)+x_move(1) x(2)+x_move(2)],[y(1)+y_move(1) y(2)+y_move(2)],'Color','k','Linewidth',2);
        line([x(2)+x_move(2) x(3)+x_move(3)],[y(2)+y_move(2) y(3)+y_move(3)],'Color','k','Linewidth',2);
        line([x(3)+x_move(3) x(4)+x_move(4)],[y(3)+y_move(3) y(4)+y_move(4)],'Color','k','Linewidth',2);
        line([x(4)+x_move(4) x(5)+x_move(5)],[y(4)+y_move(4) y(5)+y_move(5)],'Color','k','Linewidth',2);
        line([x(5)+x_move(5) x(6)+x_move(6)],[y(5)+y_move(5) y(6)+y_move(6)],'Color','k','Linewidth',2);
        line([x(6)+x_move(6) x(1)+x_move(1)],[y(6)+y_move(6) y(1)+y_move(1)],'Color','k','Linewidth',2);
      
   end
    
    
    
    
    
end
axis equal

% 
% if export.profiles==1
%     figure(10002)    
%     shading interp
%     colormap jet
%     colorbar
%     print('-djpeg',[name.directory_profiles, num2str(i_f) ,'without_mesh']);
%     display_mesh_light
%     print('-djpeg',[name.directory_profiles, num2str(i_f) ,'with_mesh']);
%     close(10002) 
%     figure(10002)    
%     display_mesh_light
%     print('-djpeg',[name.directory_profiles, num2str(i_f) ,'only_mesh']);
% 
% end
% 
% if sum(ismember(floor(elem.label/1000),[1 4 5]))~=0
%     figure(10001)
%     colormap jet
%     colorbar
% end
% 
% 
% if sum(ismember(floor(elem.label/1000),[0 2 3 4 5]))~=0
%     figure(10002)
%     colormap jet
%     colorbar
% %     figure(11002)
% %     colormap jet
% %     colorbar
% end
% 
