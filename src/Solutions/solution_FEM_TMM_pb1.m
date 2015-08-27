% solution_FEM_TMM_pb1.m
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


% Calcul de la solution analytique
        k_parallel=k_air*sin(theta);
        k_perp=sqrt(k_eq^2-k_parallel^2);
        k_z=k_air*cos(theta);
        %[R A B C D]
        x=0;
        M(1,4)=-1;
        M(1,5)=1;
        x=-dy;
        M(2,2)=-k_perp*exp(-j*k_perp*x);
        M(2,3)= k_perp*exp( j*k_perp*x);
        M(2,4)= k_z*exp(-j*k_z*x);
        M(2,5)=-k_z*exp( j*k_z*x);
        M(3,2)=K_eq_til*k_eq^2*exp(-j*k_perp*x);
        M(3,3)=K_eq_til*k_eq^2*exp( j*k_perp*x);
        M(3,4)=-K_0*k_air^2*exp(-j*k_z*x);
        M(3,5)=-K_0*k_air^2*exp( j*k_z*x);
        x=-dy-delta_y_MMT;
        FF(4,1)=-k_z*exp(-j*k_z*x);
        M(4,1)= k_z*exp( j*k_z*x);
        M(4,2)= k_perp*exp(-j*k_perp*x);
        M(4,3)=-k_perp*exp( j*k_perp*x);
        FF(5,1)=K_0*k_air^2*exp(-j*k_z*x);
        M(5,1)= K_0*k_air^2*exp( j*k_z*x);
        M(5,2)=-K_eq_til*k_eq^2*exp(-j*k_perp*x);
        M(5,3)=-K_eq_til*k_eq^2*exp( j*k_perp*x);
        
        X=M\FF;
        R_analytic(i_theta)=-X(1);
        abs_analytic(i_theta)=1-abs(R_analytic(i_theta))^2;