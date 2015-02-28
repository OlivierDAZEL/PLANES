% int_edge_vectorielle.m
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

function f=int_edge_1vectorielle(k_1,a,b,centres)

% return a matrix
% f(i)=\int_a^b e^{(k_i*(x-centres(1)))}

a_1=a-centres(1:2,1);

t=(b-a)/norm(b-a);

k1t=k_1.'*t;

temp=(k1t.');

[i_temp,j_temp]=find(abs(temp)<1e-6);
temp(i_temp,j_temp)=0;

temp0=temp.*temp==0;
tempn0=temp~=0;

temp=temp+temp0; % Pour la division

f=norm(b-a)*temp0+tempn0.*((exp(norm(b-a)*temp)-1)./temp);

%%%%
%size(f);
%size(exp((k_1.'*a_1)*ones(1,length(k_2))));
%ones(length(k_1),1);
%(k_2.'*a_2);



f=(f.').*exp((k_1.'*a_1));
