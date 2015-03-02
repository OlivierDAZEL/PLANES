% create_wave_vectors.m
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


k_air=omega/air.c;
k_x=k_air*sin(theta_inc);
k_z=k_air*cos(theta_inc);

nb_Bloch_waves=floor((period/(2*pi))*(3*real(k_air)-k_x))+5;
nb_Bloch_waves=0

if nb_R~=0
    nb_R=2*nb_Bloch_waves+1;
end
if nb_T~=0
    nb_T=2*nb_Bloch_waves+1;
end

temp=[];
temp(1:2:2*nb_Bloch_waves)=1:nb_Bloch_waves;
temp(2:2:2*nb_Bloch_waves+1)=-(1:nb_Bloch_waves);
temp=[0 temp];

vec_k_x=k_x+temp*(2*pi/period);
vec_k_x_t=k_x+temp*(2*pi/period);

vec_k_z=sqrt(k_air^2-vec_k_x.^2);
vec_k_z=real(vec_k_z)-1i*(imag(vec_k_z));

vec_k_z_t=sqrt(k_air^2-vec_k_x_t.^2);
vec_k_z_t=real(vec_k_z_t)-1i*imag(vec_k_z_t);




