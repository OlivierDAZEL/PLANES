% evaluate_polynom_2D.m
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




function f=evaluate_polynom_2D(P,x,y)


nb_v=length(x);

if nb_v==1
    nx=size(P,2);
    ny=size(P,1);
    
    Mx=(ones(ny,1)*(0:nx-1));
    My=(ones(nx,1)*(0:ny-1))';
    f=sum(sum(P.*((x.^Mx).*(y.^My))));
else
    
    nb_x=size(P,2);
    nb_y=size(P,1);
    % number of values
    PP=P(:,:,ones(nb_v,1));
    
    XX=kron(x,ones(nb_y,nb_x));
    XX=reshape(XX,[nb_y nb_x nb_v]);
    YY=kron(y,ones(nb_y,nb_x));
    YY=reshape(YY,[nb_y nb_x nb_v]);
    
    power_x=reshape(kron(ones(1,nb_v),(ones(nb_y,1)*(0:nb_x-1))),[nb_y nb_x nb_v]);
    power_y=reshape(kron(ones(1,nb_v),(ones(nb_x,1)*(0:nb_y-1))'),[nb_y nb_x nb_v]);
    
    f=PP.*(XX.^power_x).*(YY.^power_y);
    f=sum(f,1);
    f=sum(f,2);
    f=reshape(f,[1 nb_v]);
 end




