% DGM_resolution.m
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

for ii=1:nb_theta
    Shift_fluid((ii-1)*3+(1:3),ii)=1;
end


tic
for i_f=1:abs(nb_frequencies)
    DGM_progress=100*i_f/abs(nb_frequencies)
    
    freq=vec_freq(i_f);
    omega=2*pi*freq;
    k_air=omega/air.c;
    
    % Construction of the linear system
    
    A= sparse(nb_dof_DGM,nb_dof_DGM);
    F=zeros(nb_dof_DGM,1);
    
    
    for ie=1:nb_internal
        
        edge_internal
        
    end
    
    
    
    for ie=1:nb_dirichlets
        switch dirichlets(ie,4)
            case {1}
                boundary_rigid_wall
            case {5}
                boundary_sliding
            case {6}
                boundary_bonded
            case {9}
                boundary_rigid_wall_PML
        end
    end
    
    for ie=1:nb_loads
        boundary_normal_velocity
    end
    
    
    disp('Resolution of the system')
    X=A\F;
    
    %trace_DGM_y
    
    %info_DGM
    
end