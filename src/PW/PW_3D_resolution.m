% PW_resolution.m
%
% Copyright (C) 2016 < Olivier DAZEL <olivier.dazel@univ-lemans.fr> >
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

compute_number_PW_3D



for i_f=1:abs(frequency.nb)
    omega=2*pi*frequency.vec(i_f);
    k_air=omega/air.c;
    k_x=k_air*sin(data_model.theta_x);
    k_y=k_air*sin(data_model.theta_y);
    
    for i_m=1:nb_multilayers_3D
        
        Mat_PW=build_global_PW_matrices_3D(k_x,k_y,omega,multilayer_3D(:,i_m),termination_3D(i_m),nb_layers_3D(i_m),nb_amplitudes_3D(i_m),n_w_3D(:,i_m),k_air,air);
        

        F_PW=-Mat_PW(:,1);
        Mat_PW(:,1)=[];
        
        X_PW=Mat_PW\F_PW;
        
        abs_PW_3D(i_f,i_m)=1-abs(X_PW(1))^2;
        rflx_PW_3D(i_f,i_m)=X_PW(1);
        if termination~=0
            TL_PW_3D(i_f,i_m)=-20*log10(abs(X_PW(end)));
        end
        if exist([name.project_full '_postprocess'])==2
            eval([name.project_full '_postprocess'])
        end
        
    end
    
    
end
