% plot_sol_TR6_y.m
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


if sum(ismember(floor(element_label/1000),[1 4 5]))~=0
    figure(2002)
    hold on
    title('Solid displacement')
    xlabel('y')
    ylabel('abs(u)')
    figure(4002)
    hold on
    title('Solid displacement')
    xlabel('y')
    ylabel('angle(u)')
end


if sum(ismember(floor(element_label/1000),[0 2 3 4 5]))~=0
    figure(2010)
    hold on
    title('Pressure')
    xlabel('y')
    ylabel('abs(P)')
    figure(4010)
    hold on
    title('Pressure')
    xlabel('y')
    ylabel('angle(P)')    
end



for ie=1:nb.elements
    
    
    if ismember(floor(element_label(ie)/1000),[1 4 5])
        c=sqrt(sol(3*(elements(ie,:)-1)+1).^2+sol(3*(elements(ie,:)-1)+2).^2);
        y=nodes(elements(ie,:),2);
        figure(2010)
        plot(y,abs(c),'r.');
        figure(4010)
        plot(y,angle(c),'r.');        
    end
    
    if ismember(floor(element_label(ie)/1000),[0 2 4 5 8])
        c=sol(3*(elements(ie,:)-1)+3);
        y=nodes(elements(ie,:),2);
        figure(2010)
        plot(y,abs(c),'r.');    
        figure(4010)
        plot(y,angle(c),'r.');         
    end
    
end


if export_profiles==1
    shading interp
    print('-djpeg',[name_directory_profiles, num2str(i_f)]);
end

% if sum(ismember(floor(element_label/1000),[1 4 5]))~=0
%     close(10001)
% end
% 
% 
% if sum(ismember(floor(element_label/1000),[0 2 3 4 5]))~=0
%     close(10002)
% end

