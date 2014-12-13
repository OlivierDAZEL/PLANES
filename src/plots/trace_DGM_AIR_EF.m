% trace_DGM_AIR_EF.m
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




for ie=1:nb_elements
    coord_elem=vcor(kconec(ie,:),1:2)';
    x_centre=centre_element(ie);

    e_edge=ie;
    parameter_element
    
        q=X_DGM_fortran(dof_start_element(ie)+[0:ondes_element(ie)*nb_theta-1]);
    Phi_elem=zeros(3,3);
    for i_thetaphi=1:nb_theta
        theta_phi=vec_theta(i_thetaphi);
        n_phi=[cos(theta_phi);sin(theta_phi)];
        Phi_e=Phi_fluid(cos(theta_phi),sin(theta_phi),Z_e);
        Phi_elem(:,1)=Phi_elem(:,1)+Phi_e*exp(-j*k_e*(n_phi'*(coord_elem(:,1)-x_centre)))*q(1+ondes_element(ie)*(i_thetaphi-1));
        Phi_elem(:,2)=Phi_elem(:,2)+Phi_e*exp(-j*k_e*(n_phi'*(coord_elem(:,2)-x_centre)))*q(1+ondes_element(ie)*(i_thetaphi-1));
        Phi_elem(:,3)=Phi_elem(:,3)+Phi_e*exp(-j*k_e*(n_phi'*(coord_elem(:,3)-x_centre)))*q(1+ondes_element(ie)*(i_thetaphi-1));
    end

    x=[coord_elem(1,:)];
    y=[coord_elem(2,:)];
    c=[Phi_elem(1,1);Phi_elem(1,2);Phi_elem(1,3)];
    
    c_trace=(c);

    figure(1111)
    hold on
    patch(x,y,abs((c)));
  %  colorbar


        % Trac? des figures pour le cas 1D
        figure(1000)
        for i_fig=1
            %subplot(2,2,i_fig)
            hold on
            plot(x,abs([Phi_elem(i_fig,1);Phi_elem(i_fig,2);Phi_elem(i_fig,3)]),'.')
        end
%         figure(1001)
%         for i_fig=1
%             %subplot(2,2,i_fig)
%             hold on
%             plot(x,angle([Phi_elem(i_fig,1);Phi_elem(i_fig,2);Phi_elem(i_fig,3)]),'.')
%         end
        
        
        
%         figure(1001)
%         for i_fig=1:3
%             subplot(2,2,i_fig)
%             hold on
%             plot(x,angle([Phi_elem(i_fig,1);Phi_elem(i_fig,2);Phi_elem(i_fig,3)]),'.')
%         end



end