% create_vector_waves.m
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








for ie=1:nb_elements
    if floor(element_label(ie)/1000)==0
        ondes_element(ie)=1;   % Air
    elseif floor(element_label(ie)/1000)==1
        ondes_element(ie)=2;   % Milieu solide elastique
    elseif floor(element_label(ie)/1000)==2
        ondes_element(ie)=1;   % Fluide equivalent
     elseif floor(element_label(ie)/1000)==3
        ondes_element(ie)=1;   % Limp
      elseif floor(element_label(ie)/1000)==4
        ondes_element(ie)=3;   % Fluide equivalent
      elseif floor(element_label(ie)/1000)==8
        ondes_element(ie)=1;   % Fluide equivalent  
    else
        disp('Subroutine import_mesh')
        disp('Unknwon type of element')
        break
    end
end

angle_element=nb_theta*ones(nb_elements,1);

dof_start_element(1)=1; %
for ie=2:nb_elements
   dof_start_element(ie)=dof_start_element(ie-1)+ondes_element(ie-1)*angle_element(ie-1);
end

