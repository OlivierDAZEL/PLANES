% check_dirichlet.m
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


isvalidof=zeros(3*nb_nodes,1);

for ie=1:nb_elements

	typ=floor(element_label(ie)/1000);
	switch typ
		case {0,2,3,8} %! Acoustic/EF/limep
					isvalidof(3*(elements(ie,1:6)-1)+3)=1;
        case {1} %! Elastic solid
					isvalidof(3*(elements(ie,1:6)-1)+1)=1;
					isvalidof(3*(elements(ie,1:6)-1)+2)=1;

         case {4,5}	%! PEM
					isvalidof(3*(elements(ie,1:6)-1)+1)=1;
					isvalidof(3*(elements(ie,1:6)-1)+2)=1;
					isvalidof(3*(elements(ie,1:6)-1)+3)=1;
        otherwise
				disp('Attention element sans nature connue')
                stop
    end
end   



for ie=1:nb_dirichlets
	if (dirichlets(ie,4)==5) % Sliding
		xx=abs(nodes(dirichlets(ie,1),1)-nodes(dirichlets(ie,2),1));
		yy=abs(nodes(dirichlets(ie,1),2)-nodes(dirichlets(ie,2),2));
		if (xx>yy) 
			isvalidof(3*(dirichlets(ie,1)-1)+2)=0;
			isvalidof(3*(dirichlets(ie,2)-1)+2)=0;
			isvalidof(3*(dirichlets(ie,6)-1)+2)=0;
        else
			isvalidof(3*(dirichlets(ie,1)-1)+1)=0;
			isvalidof(3*(dirichlets(ie,2)-1)+1)=0;
			isvalidof(3*(dirichlets(ie,6)-1)+1)=0;
        end
    end
	if (dirichlets(ie,4)==6) % Bonded
        isvalidof(3*(dirichlets(ie,1)-1)+1)=0;
		isvalidof(3*(dirichlets(ie,1)-1)+2)=0;
		isvalidof(3*(dirichlets(ie,2)-1)+1)=0;
		isvalidof(3*(dirichlets(ie,2)-1)+2)=0;
		isvalidof(3*(dirichlets(ie,6)-1)+1)=0;
		isvalidof(3*(dirichlets(ie,6)-1)+2)=0;
    end
	if (dirichlets(ie,4)==1) 
		isvalidof(3*(dirichlets(ie,1)-1)+1)=0;
		isvalidof(3*(dirichlets(ie,1)-1)+2)=0;
		isvalidof(3*(dirichlets(ie,2)-1)+1)=0;
		isvalidof(3*(dirichlets(ie,2)-1)+2)=0;
		isvalidof(3*(dirichlets(ie,6)-1)+1)=0;
		isvalidof(3*(dirichlets(ie,6)-1)+2)=0;
    end	
end





dof_A=zeros(3*nb_nodes,1);


dof_back=[];



itemp=0;
for ii=1:3*nb_nodes
	if (isvalidof(ii)) 
		itemp=itemp+1;
		dof_A(ii)=itemp;
		dof_back(itemp)=ii;
    end
end

nb_dof_FEM=itemp;
list_dof_valid=find(isvalidof);

