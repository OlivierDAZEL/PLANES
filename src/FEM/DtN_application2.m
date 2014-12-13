% DtN_application2.m
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


if (nb_R+nb_T)>0
        A_2(nb_dof_FEM+nb_R+nb_T,nb_dof_FEM+nb_R+nb_T)=0;
        F_2(nb_dof_FEM+nb_R+nb_T,1)=0;
end
    
    DtN_porous_R=0;
    DtN_porous_T=0;
    DtN_elas_R=0;
    DtN_elas_T=0;
    
    for ie=1:nb_loads  
        typ=loads(ie,4);
        switch typ
            case {10}
                DtN_porous_R=1;
                excitation_10000
            case {11}
                DtN_elas_R=1;
                excitation_20000
            case {21}
                DtN_elas_T=1;
                excitation_21000
            case {12}
                DtN_porous_R=1;
                excitation_50000
            otherwise
                disp('Unknown load')
                stop
        end
    end
    
 
    
    if (nb_media.eqf*DtN_porous_R)~=0
        %    DtN Biot 1998
        if (nb_R~=0)
            F_2(nb_dof_FEM+nb_Bloch_waves+1)=period;
            for i_R=1:nb_R
                 A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)=A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)-period;
            end
        end
    end
   
        if (nb_media.limp*DtN_porous_R)~=0
        %    DtN Biot 1998
        if (nb_R~=0)
            F_2(nb_dof_FEM+nb_Bloch_waves+1)=period;
            for i_R=1:nb_R
                 A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)=A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)-period;
            end
        end
    end
    
    
    
    
if (nb_media.pem98*DtN_porous_R)~=0
        %    DtN Biot 1998
        if (nb_R~=0)
            F_2(nb_dof_FEM+nb_Bloch_waves+1)=period;
            for i_R=1:nb_R
                 A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)=A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)-period;
            end
        end
    end  
    
    
    
    if (nb_media.pem01*DtN_porous_R)~=0
        %    DtN Biot 2001
        if (nb_R~=0)
            F_2(nb_dof_FEM+1)=period;
            for i_R=1:nb_R
                A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)=A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)-period;
            end
        end
    end
    
    
    
    
    
    if (nb_media.elas*DtN_elas_R)~=0
        %    DtN Biot elas
        if (nb_R~=0)
            F_2(nb_dof_FEM+nb_Bloch_waves+1)=F_2(nb_dof_FEM+nb_Bloch_waves+1)+period*(1i*k_z)/(rho_0*omega^2);
            for i_R=1:nb_R
                A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)=A_2(nb_dof_FEM+i_R,nb_dof_FEM+i_R)-...
                    period*(1i*vec_k_z(i_R))/(rho_0*omega^2);
            end
        end
        if (nb_T~=0)
            for i_T=1:nb_T
                A_2(nb_dof_FEM+nb_R+i_T,nb_dof_FEM+nb_R+i_T)=A_2(nb_dof_FEM+nb_R+i_T,nb_dof_FEM+nb_R+i_T)+...
                    period*(1i*vec_k_z(i_T))/(rho_0*omega^2);
            end
        end
    end