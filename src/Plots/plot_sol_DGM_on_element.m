% plot_sol_DGM_on_element.m
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


coord_elem=nodes(nonzeros(elem.nodes(ie,:)),1:2)';
x_centre=mean(nodes(nonzeros(elem.nodes(ie,:)),1:2))';

q=X(dof_start_element(ie)+[0:ondes_element(ie)*theta_DGM.nb-1]);


e_edge=ie;
parameter_element

if elem.model(ie)==10
    Phi_elem=zeros(3,3);
else
    Phi_elem=zeros(3,4);
end
for i_thetaphi=1:theta_DGM.nb
    theta_phi=vec_theta(i_thetaphi);
    n_phi=[cos(theta_phi)*tau_x;sin(theta_phi)*tau_y];
    Phi_e=Phi_fluid(cos(theta_phi),sin(theta_phi),Z_e);
    Phi_elem(:,1)=Phi_elem(:,1)+Phi_e*exp(-1j*k_e*(n_phi.'*(coord_elem(:,1)-x_centre)))*q(1+ondes_element(ie)*(i_thetaphi-1));
    Phi_elem(:,2)=Phi_elem(:,2)+Phi_e*exp(-1j*k_e*(n_phi.'*(coord_elem(:,2)-x_centre)))*q(1+ondes_element(ie)*(i_thetaphi-1));
    Phi_elem(:,3)=Phi_elem(:,3)+Phi_e*exp(-1j*k_e*(n_phi.'*(coord_elem(:,3)-x_centre)))*q(1+ondes_element(ie)*(i_thetaphi-1));
    if elem.model(ie)==11
        Phi_elem(:,4)=Phi_elem(:,4)+Phi_e*exp(-1j*k_e*(n_phi.'*(coord_elem(:,4)-x_centre)))*q(1+ondes_element(ie)*(i_thetaphi-1));
    end
end

y=[coord_elem(2,:)];
c=[Phi_elem(3,:)];
figure(2010)
plot(y,abs(c),'b.');
figure(4010)
plot(y,angle(c),'b.');







