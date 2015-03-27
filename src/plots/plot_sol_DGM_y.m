% trace_DGM_y.m
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
%                 https://github.com/OlivierDAZEL/PLANES
% or find more details on Olivier's webpage
%                  http://perso.univ-lemans.fr/~odazel/
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
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


figure(10002)
hold on
% figure(10003)
% hold on



for ie=1:nb_elements
    coord_elem=nodes(elements(ie,:),1:2)';
    x_centre=mean(nodes(elements(ie,:),1:2))';
    
    q=X(dof_start_element(ie)+[0:ondes_element(ie)*nb_theta-1]);
    
    switch floor(element_label(ie)/1000)
        case {0,2,8}
            e_edge=ie;
            parameter_element
            
            Phi_elem=zeros(3,3);
            for i_thetaphi=1:nb_theta
                theta_phi=vec_theta(i_thetaphi);
                n_phi=[cos(theta_phi)*tau_x;sin(theta_phi)*tau_y];
                Phi_e=Phi_fluid(cos(theta_phi),sin(theta_phi),Z_e);
                Phi_elem(:,1)=Phi_elem(:,1)+Phi_e*exp(-1j*k_e*(n_phi.'*(coord_elem(:,1)-x_centre)))*q(1+ondes_element(ie)*(i_thetaphi-1));
                Phi_elem(:,2)=Phi_elem(:,2)+Phi_e*exp(-1j*k_e*(n_phi.'*(coord_elem(:,2)-x_centre)))*q(1+ondes_element(ie)*(i_thetaphi-1));
                Phi_elem(:,3)=Phi_elem(:,3)+Phi_e*exp(-1j*k_e*(n_phi.'*(coord_elem(:,3)-x_centre)))*q(1+ondes_element(ie)*(i_thetaphi-1));
            end
            
            y=[coord_elem(2,:)];
            
            figure(10002)
            hold on
            c=[Phi_elem(3,1);Phi_elem(3,2);Phi_elem(3,3)];
            plot(y,abs(c),'r+')
%             figure(10003)
%             hold on
%             c=[Phi_elem(3,1);Phi_elem(3,2);Phi_elem(3,3)];
%             plot(y,angle(c),'r+')
            
            
            
        case {4,5}
            
            k_e=[delta_1 delta_2 delta_3];
            Phi_elem=zeros(8,3);
            for i_thetaphi=1:nb_theta
                theta_phi=vec_theta(i_thetaphi);
                n_phi=[cos(theta_phi);sin(theta_phi)];
                Phi_e=Phi_Biot(cos(theta_phi),sin(theta_phi),delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
                for i_onde=1:3
                    Phi_elem(:,1)=Phi_elem(:,1)+Phi_e(:,i_onde)*exp(-1j*k_e(i_onde)*(n_phi'*(coord_elem(:,1)-x_centre)))*q(i_onde+ondes_element(ie)*(i_thetaphi-1));
                    Phi_elem(:,2)=Phi_elem(:,2)+Phi_e(:,i_onde)*exp(-1j*k_e(i_onde)*(n_phi'*(coord_elem(:,2)-x_centre)))*q(i_onde+ondes_element(ie)*(i_thetaphi-1));
                    Phi_elem(:,3)=Phi_elem(:,3)+Phi_e(:,i_onde)*exp(-1j*k_e(i_onde)*(n_phi'*(coord_elem(:,3)-x_centre)))*q(i_onde+ondes_element(ie)*(i_thetaphi-1));
                end
            end
            
            
            y=[coord_elem(2,:)];
            
            figure(10002)
            hold on
            c=[Phi_elem(8,1);Phi_elem(8,2);Phi_elem(8,3)];
            plot(y,abs(c),'r+')
            
            figure(10001)
            hold on
            c=[Phi_elem(2,1);Phi_elem(2,2);Phi_elem(2,3)]/(j*omega);
            plot(y,abs(c),'r+')
            
    end
    
    
    
end