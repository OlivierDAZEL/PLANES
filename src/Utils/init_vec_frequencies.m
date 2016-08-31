% init_vec_frequencies.m
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

if length(frequency.nb)==1
    if(frequency.nb>0)
        if frequency.nb==1
            frequency.vec=frequency.min;
        else
            frequency.vec=linspace(frequency.min,frequency.max,frequency.nb);
        end
    else
        frequency.vec=logspace(log10(frequency.min),log10(frequency.max),abs(frequency.nb));
    end
else % Case of complex frequency
    temp_1=linspace(frequency.min,frequency.max,frequency.nb(1));
    temp_2=linspace(frequency.min_imag,frequency.max_imag,frequency.nb(2));
    frequency.vec=[];
    for ii=1:frequency.nb(2)
        frequency.vec=[frequency.vec temp_1+1j*temp_2(ii)];
    end
    frequency.nb=frequency.nb(1)*frequency.nb(2);
end