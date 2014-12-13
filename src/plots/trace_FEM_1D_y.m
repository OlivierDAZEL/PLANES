% trace_FEM_1D_y.m
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



x=[];
y=[];
% figure(567)
% hold on
% for ie=1:nb_elements
% 
%     c=X_EF_fortran(3*kconec(ie,:));
%     vertices=[vcor(kconec(ie,:),:)'];    
%     plot(vertices(2,:),abs(c),'.');
%      
% end
% title('Pression EF (module)')


% hold on
% xx=linspace(0,1,200);
% plot(xx,2*abs(cos(omega*(xx-1)/c_0)),'r');





figure(5679)
hold on
for ie=1:nb_elements

    c=X_EF_fortran(3*elements(ie,1:6));
    vertices=[nodes(elements(ie,1:6),:)'];    
    plot(vertices(2,:),abs(c),'.');
     
end
title('Pression EF (module)')




figure(5678)
hold on
for ie=1:nb_elements

    c=X_EF_fortran(3*(elements(ie,1:6)-1)+2);
    vertices=[nodes(elements(ie,1:6),:)'];    
    plot(vertices(2,:),abs(c),'.');
     
end
title('u_y EF (module)')












% figure(56780)
% hold on
% for ie=1:nb_elements
% 
%     c=X_EF_fortran(3*(kconec(ie,:)-1)+2);
%     vertices=[vcor(kconec(ie,:),:)'];    
%     plot(vertices(2,:),angle(c),'.');
%      
% end
% title('u_y EF (phase)')
% 
% 
% figure(56790)
% hold on
% for ie=1:nb_elements
% 
%     c=X_EF_fortran(3*kconec(ie,:));
%     vertices=[vcor(kconec(ie,:),:)'];    
%     plot(vertices(2,:),angle(c),'.');
%      
% end
% title('Pression EF (phase)')
