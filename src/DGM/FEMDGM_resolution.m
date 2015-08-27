% FEMDGM_resolution.m
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

for ii=1:nb_thetaDGM
    Shift_fluid((ii-1)*3+(1:3),ii)=1;
    % Biot wave 1
    Shift_Biot((ii-1)*8*3+(1:8),1+(ii-1)*3)=1;
    % Biot wave 2
    Shift_Biot((ii-1)*8*3+8+(1:8),2+(ii-1)*3)=1;
    % Biot wave 3
    Shift_Biot((ii-1)*8*3+16+(1:8),3+(ii-1)*3)=1;
end

tic
for i_f=1:abs(nb_frequencies)
    DGM_progress=100*i_f/abs(nb_frequencies)
    
    freq=vec_freq(i_f);
    omega=2*pi*freq;
    k_air=omega/air.c;
    if nb.R~=0
        create_wave_vectors
        delta_Bloch=exp(-1i*k_x*period);
    end
    Mat_parameter=initialize_Mat_parameter(index_label,index_element,air,omega);
    
    % Construction of the linear system
    
    A= sparse(nb.dof_FEM+nb.dof_DGM+nb.R,nb.dof_FEM+nb.dof_DGM+nb.R);
    F=zeros(nb.dof_FEM+nb.dof_DGM+nb.R,1);
    
    for ie=1:nb.internal
            edge_internal
    end
    
    
    for ie=1:nb.dirichlets
        if element_model(dirichlets(ie,3))
            switch dirichlets(ie,4)
                case {1}
                    %disp('Lancement boundary_rigid_wall')
                    boundary_rigid_wall
                case {5}
                    boundary_sliding
                case {6}
                    boundary_bonded
                case {9}
                    boundary_rigid_wall_PML
            end
        end
    end
    
    
    
    for ie=1:nb.loads
        if element_model(loads(ie,3))==1
            %disp('Lancement boundary_normal_velocity')
            switch loads(ie,4)
                case {3}
                    boundary_normal_displacement
                case{10}
                    boundary_10
            end
        end
    end
    
    
    for ie=1:nb.periodicity
        if element_model(loads(ie,3))==1
            edge_periodicity
        end
    end
    
   if (nb.media.acou~=0)
        A(1:nb.dof_FEM,1:nb.dof_FEM)=A(1:nb.dof_FEM,1:nb.dof_FEM)+(H_acou/(air.rho*omega^2)-Q_acou/(air.K));
   end
    
        
    disp('Resolution of the system')
    X=A\F;
    sol=[];
    sol(dof_back)=X(1:nb.dof_FEM);
   
    
    if profiles.on==1
        disp('plotting the solution')
        if profiles.y==1
            plot_sol_DGM_y_Q
            plot_sol_H12_y
        end
        if profiles.custom~=0
            eval(['plot_sol_DGM_custom_' , num2str(profiles.custom)]);
        end
    end
    
    if nb.R~=0
        reflex_DGM(i_f)=X(end);
    end
    %info_DGM
    
end