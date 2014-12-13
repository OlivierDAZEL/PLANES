% periodicity_condition_application.m
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


%disp('Applying periodicity condtions')



edge_left= find(periodicity(:,4)==98);
edge_right=find(periodicity(:,4)==99);

node_left=unique([periodicity(edge_left,1);periodicity(edge_left,2);periodicity(edge_left,6)]);
[temp,i_left]=sort(nodes(node_left,2));
node_left=node_left(i_left);

node_right=unique([periodicity(edge_right,1);periodicity(edge_right,2);periodicity(edge_right,6)]);
[temp,i_right]=sort(nodes(node_right,2));
node_right=node_right(i_right);

delta=exp(-1i*k_x*period);

dof_left= dof_A([3*(node_left -1)+1;3*(node_left -1)+2;3*(node_left -1)+3]);
dof_right=dof_A([3*(node_right-1)+1;3*(node_right-1)+2;3*(node_right-1)+3]);

dof_left= dof_left (find(dof_left));
dof_right=dof_right(find(dof_right));


for ii=1:length(dof_left)
    
        A(:,dof_left(ii))=A(:,dof_left(ii))+delta*A(:,dof_right(ii));
        A(:,dof_right(ii))=0;

        
        A(dof_left(ii),:)=A(dof_left(ii),:)+A(dof_right(ii),:)/delta;
        F(dof_left(ii))=F(dof_left(ii))+F(dof_right(ii))/delta;

        A(dof_right(ii),:)=0;
        
        F(dof_right(ii))=0;
        A(dof_right(ii),dof_left(ii))=delta;
        A(dof_right(ii),dof_right(ii))=-1;
        
    
end