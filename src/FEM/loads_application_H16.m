% loads_application_H16.m
%
% Copyright (C) 2015 < Olivier DAZEL <olivier.dazel@univ-lemans.fr> >
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
    
    
    
    node_1=loads(ie,1);
    node_2=loads(ie,2);
    
    x1=nodes(loads(ie,1),1);
    y1=nodes(loads(ie,1),2);
    x2=nodes(loads(ie,2),1);
    y2=nodes(loads(ie,2),2);
    length_edge=sqrt((x2-x1)^2+(y2-y1)^2);
    lx=length_edge;
    
    load_Hermite_2D

    index_p=dof_A(p_H(elements(loads(ie,3),:)));
    indice_test=index_p([1 2 6 5]);
    indice_champs_0=index_p([1 2 6 5]);
    indice_champs_1=index_p([3 4 16 7]);
    
    switch typ
        case {3}
            if sort(loads(ie,1:2))==sort(elements(loads(ie,3),1:2))
                y=0;

                for i_test=1:4
                    eval(['Psi_test=Psi_',num2str(i_test),'_x'])
                    F(indice_test(i_test))=integrate_polynom(Psi_test,lx);                    
                end
    
            else    
               aezezaezaezeazezaeazezaeazezaeza 
                
                
            end
            
    end
end


