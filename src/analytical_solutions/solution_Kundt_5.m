% solution_Kundt.m
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


x_air=linspace(-model_data.ly,0,1000);
x_trace=x_air+model_data.ly;
%u(x_air)=A cos(k x_air)+B sin (k x_air)
%u(0)=0 -> A=0
%v(-ly)=1 -> 
B=1/(-1j*omega*sin(k_air*model_data.ly));
u_air=B*sin(k_air*x_air);
v_air=1j*omega*u_air;
% p_air=-K_air*u'(x_air)
p_air=-B*(air.K)*k_air*cos(k_air*x_air);


if profiles.on~=0
    

    if profiles.y~=0
        figure(2002)
        hold on
        plot(x_air+model_data.ly,abs(v_air),'b')
        figure(4002)
        hold on
        plot(x_air+model_data.ly,angle(v_air),'b')
        figure(2010)
        hold on
        plot(x_trace,abs(p_air),'b')
        figure(4010)
        hold on
        plot(x_air+model_data.ly,angle(p_air),'b')
    end
    

end
