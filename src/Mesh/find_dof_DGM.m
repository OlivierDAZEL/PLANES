% find_dof_DGM.m
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


temp=find(ismember(elem.model,[10 11]));


if length(temp)>0
    for ie=1:length(temp)
        switch floor(elem.label(ie)/1000)
            case {0,2,3,8}
                ondes_element(temp(ie))=1;
            case {1}
                ondes_element(temp(ie))=2;
            case {4,5}
                ondes_element(temp(ie))=3;
            otherwise
                pas_de_nature_connue
        end
    end
    dof_start_element(temp(1))=nb.dof_FEM+1;
    nb.dof_DGM=dof_start_element(temp(1))+data_model.theta_DGM.nb*ondes_element(temp(1))-1-nb.dof_FEM;
    
    if length(temp)>1
        for ie=2:length(temp)
            dof_start_element(temp(ie))=dof_start_element(temp(ie-1))+data_model.theta_DGM.nb*ondes_element(temp(ie-1));
        end
        nb.dof_DGM=dof_start_element(temp(ie))+data_model.theta_DGM.nb*ondes_element(temp(ie))-1-nb.dof_FEM;
    end
    
    nb.dof_radiative=0;
    for ie=1:nb.radiative
        
    end
    
    
    
    
    
    
        
    for ii=1:data_model.theta_DGM.nb
        Shift_fluid((ii-1)*3+(1:3),ii)=1;
        % Biot wave 1
        Shift_Biot((ii-1)*8*3+(1:8),1+(ii-1)*3)=1;
        % Biot wave 2
        Shift_Biot((ii-1)*8*3+8+(1:8),2+(ii-1)*3)=1;
        % Biot wave 3
        Shift_Biot((ii-1)*8*3+16+(1:8),3+(ii-1)*3)=1;
    end
else
    nb.dof_DGM=0;
end

