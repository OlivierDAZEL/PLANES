% compute_number_PW_TMM.m
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

for i_m=1:nb_multilayers
    for ii=1:multilayer(1,i_m).nb
        floor(multilayer(ii,i_m).mat/1000)
        switch floor(multilayer(ii,i_m).mat/1000)
            case 1
                n_w(ii,i_m)=4;
            case {0 2 3}
                n_w(ii,i_m)=2;
            case {4 5}
                n_w(ii,i_m)=6;
        end
    end
    nb_amplitudes(i_m)=sum(n_w(:,i_m));
    if multilayer(1,i_m).termination~=0
        nb_amplitudes(i_m)=nb_amplitudes(i_m)+1;
    end
end
