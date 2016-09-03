% properties_jca.m
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


switch porous_model.eqf
    case{'JCA'}
        properties_JCA
    case{'JCA_aniso'}
        properties_JCA_aniso
    case{'JKD_Lafarge'}    
    case{'Pride_Lafarge'}
    case{'Delany_Bazley'}
    case{'Horoshenkov'}
    case{'Wilson'}  
    case{'air'}    
    case{'custom'} 
        rho_eq_til=air.rho;
        K_eq_til=air.K*0.5^2;
        %K_eq_til=air.K;
    otherwise
        PLANES_INTERRUPTED_in_properties_eqf
end

if typ_mat==3
    eqf2limp
end




        
    