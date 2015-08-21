% plot_sol_H12_on_element.m
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



vertices=[nodes(nonzeros(elem.nodes(ie,:)),:)'];
%changement de base (cart to triangle)


if ismember(floor(elem.label(ie)/1000),[0 4 5 8])
    figure(10002)
    c=sol(p_H12(nonzeros(elem.nodes(ie,:))));

    
    
    for i_faces=1:size(faces,2)
        figure(10002)
        patch(vertices(1,:),vertices(2,:),abs(c(1:3:end)));
    end
        for i_faces=1:size(faces,2)
        figure(11002)
        patch(vertices(1,:),vertices(2,:),angle(c(1:3:end)));
    end
    
end







