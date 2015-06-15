% loads_application_H16.m
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

for ie=1:nb.loads
    typ=floor(loads(ie,4));
    
    x1=nodes(loads(ie,1),1);
    y1=nodes(loads(ie,1),2);
    x2=nodes(loads(ie,2),1);
    y2=nodes(loads(ie,2),2);
    length_edge=sqrt((x2-x1)^2+(y2-y1)^2);
    
    
    xx=abs(x2-x1);
	yy=abs(y2-y1);

    
    
    
    switch typ
        case {3}
            index_force=p_H(loads(ie,1:2));
            index_force=index_force(1:4:end);
            index_force=dof_A(index_force);        
            index_F_elem=find(index_force);
            index_F_global=index_force(index_F_elem); 
            F(index_F_global)=F(index_F_global)-length_edge/2;
            
            index_force=p_H(loads(ie,1:2));
            if xx>yy
                index_force=index_force(2:4:end);
            else 
                index_force=index_force(3:4:end);
            end
            index_force=dof_A(index_force);        
            index_F_elem=find(index_force);
            index_F_global=index_force(index_F_elem); 
            F(index_F_global(1))=F(index_F_global(1))-length_edge/12;
             F(index_F_global(2))=F(index_F_global(2))+length_edge/12;
         
            
    end
end


