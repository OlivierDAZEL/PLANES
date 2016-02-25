% build_global_PW_matrices.m
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


function Mat_PW=build_global_PW_matrices(k_x,omega,multilayer,nb_amplitudes,n_w,k_air,air)


% Initialization of the matrix
Mat_PW=zeros(nb_amplitudes-1,nb_amplitudes);

% Creation of the equation
number_of_eq=0;
% Space shift for the first layer being x=0
x_interface=-multilayer(1).d;
% Loop on the layers
for i_interface=1:multilayer(1,1).nb-1
    
    %Type of media on both sides and attribution of the dof
    medium_1=multilayer(i_interface).mat;
    dof_medium_1=sum(n_w(1:i_interface-1))+(1:n_w(i_interface));
    medium_2=multilayer(i_interface+1).mat;
    dof_medium_2=sum(n_w(1:i_interface))+(1:n_w(i_interface+1));
    
    if ismember(floor(medium_1/1000),[0 2 3])
        switch floor(medium_2/1000)
            case {0 2 3}
                interface_fluid_fluid_2D
            case 1
                interface_fluid_elas_2D
            case {4 5}
                interface_fluid_PEM_2D
        end
    elseif floor(medium_1/1000)==1
        switch floor(medium_2/1000)
            case {0 2 3}
                interface_elas_fluid_2D
            case 1
                interface_elas_elas_2D
            case {4 5}
                interface_elas_PEM_2D
        end
    elseif ismember(floor(medium_1/1000),[4 5])
        switch floor(medium_2/1000)
            case {0 2 3}
                interface_PEM_fluid_2D
            case 1
                interface_PEM_elas_2D
            case {4 5}
                interface_PEM_PEM_2D
        end
    end
end
% Last interface
% the last medium is #2 of the end of the loop


if multilayer(1,1).termination==0 % Rigid backing
    switch floor(medium_2/1000)
        case {0 2 3}
            termination_rigid_fluid_2D
        case 1
            termination_rigid_elas_2D
        case {4 5}
            termination_rigid_PEM_2D
    end
else % transmission problem
    % Medium # 1= air last medium is #2
    k_z_1=sqrt(k_air^2-k_x^2);
    SV_1=State_fluid_2D(k_x,k_z_1,air.K);
    switch floor(medium_2/1000)
        case {0 2 3}
            termination_trans_fluid_2D
        case 1
            termination_trans_elas_2D
        case {4 5}
            termination_trans_PEM_2D
    end
    
end

end