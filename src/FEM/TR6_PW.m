% TR6_PW.m
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


function F = TR6_PW(longueur,k,a)

n_pg=4;


points_gauss(1)= 0.339981043584856;
points_gauss(2)=-0.339981043584856;
points_gauss(3)= 0.861136311594053;
points_gauss(4)=-0.861136311594053;

weight_gauss(1)=0.652145154862546;
weight_gauss(2)=0.652145154862546;
weight_gauss(3)=0.347854845137454;
weight_gauss(4)=0.347854845137454;


F=zeros(3,1);

for i_pg=1:n_pg
    
    xi=points_gauss(i_pg);
    
    F(1)=F(1)+(-(1-xi)*xi/2)  *exp(-1i*k*(longueur*xi/2))*(weight_gauss(i_pg));
    F(2)=F(2)+( (1+xi)*xi/2)  *exp(-1i*k*(longueur*xi/2))*(weight_gauss(i_pg));
    F(3)=F(3)+( (1-xi)*(1+xi))*exp(-1i*k*(longueur*xi/2))*(weight_gauss(i_pg));
end

F=F*(longueur/2)*exp(-1i*k*(a+longueur/2));


end