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

if sum(ismember(floor(element_label/1000),[1 4 5]))~=0
    figure(2005)
    hold on
    title('Solid displacement')
    xlabel('y')
    ylabel('abs(u)')
    figure(4005)
    hold on
    title('Solid displacement')
    xlabel('y')
    ylabel('angle(u)')
end


if sum(ismember(floor(element_label/1000),[0 2 3 4 5]))~=0
    figure(2002)
    hold on
    title('Fluid displacement')
    xlabel('y')
    ylabel('abs(u)')
    figure(4002)
    hold on
    title('Fluid displacement')
    xlabel('y')
    ylabel('angle(u)')
    figure(2010)
    hold on
    title('Pressure')
    xlabel('y')
    ylabel('abs(P)')
    figure(4010)
    hold on
    title('Pressure')
    xlabel('y')
    ylabel('angle(P)')
end



for ie=1:nb.elements
    if element_model(ie)==1
        nb_nodes_elements=length(elements(ie,:));
        vertices=[nodes(elements(ie,:),1:2)';zeros(1,nb_nodes_elements)];
        if (length(elements(ie,:)==3))
            faces = [1 2 3]';
            for ii=1:1
                [vertices, faces]=linearSubdivision(vertices, faces);
            end
       end
       x_centre=centre_element(ie,nodes,elements);
        
        e_edge=ie;
        
        q=X(dof_start_element(ie)+[0:ondes_element(ie)*nb_thetaDGM-1]);
        vec_normales=[cos(vec_theta);sin(vec_theta)]';
        
        switch floor(element_label(ie)/1000)
            case {0,2,8}
                e_edge=ie;
                parameter_element
                
                for i_point=1:size(vertices,2);
                    x_point=diag([tau_x tau_y])*vertices(1:2,i_point);
                    Delta=diag(exp(-1j*k_e*vec_normales*(x_point-x_centre)));
                    Phi_point(1:3,i_point)=sum(Phi_fluid(cos(vec_theta),sin(vec_theta),Z_e)*Delta*q,2);
                end
                figure(2010)
                for i_faces=1:size(faces,2)
                    plot(vertices(2,faces(:,i_faces)),abs(Phi_point(3,faces(:,i_faces))'),'m.');
                end
                figure(4010)
                for i_faces=1:size(faces,2)
                    plot(vertices(2,faces(:,i_faces)),angle(Phi_point(3,faces(:,i_faces))'),'m.');
                end
                figure(2002)
                for i_faces=1:size(faces,2)
                    plot(vertices(2,faces(:,i_faces)),abs(Phi_point(2,faces(:,i_faces))'),'m.');
                end
                figure(4002)
                for i_faces=1:size(faces,2)
                    plot(vertices(2,faces(:,i_faces)),angle(Phi_point(2,faces(:,i_faces))'),'m.');
                end
                
            case {4,5}
                
                k_e=[delta_1 delta_2 delta_3];
                Phi_elem=zeros(8,4);
                for i_thetaphi=1:nb_theta
                    theta_phi=vec_theta(i_thetaphi);
                    n_phi=[cos(theta_phi);sin(theta_phi)];
                    Phi_e=Phi_Biot(cos(theta_phi),sin(theta_phi),delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
                    for i_onde=1:3
                        Phi_elem(:,1)=Phi_elem(:,1)+Phi_e(:,i_onde)*exp(-1j*k_e(i_onde)*(n_phi'*(coord_elem(:,1)-x_centre)))*q(i_onde+ondes_element(ie)*(i_thetaphi-1));
                        Phi_elem(:,2)=Phi_elem(:,2)+Phi_e(:,i_onde)*exp(-1j*k_e(i_onde)*(n_phi'*(coord_elem(:,2)-x_centre)))*q(i_onde+ondes_element(ie)*(i_thetaphi-1));
                        Phi_elem(:,3)=Phi_elem(:,3)+Phi_e(:,i_onde)*exp(-1j*k_e(i_onde)*(n_phi'*(coord_elem(:,3)-x_centre)))*q(i_onde+ondes_element(ie)*(i_thetaphi-1));
                        Phi_elem(:,4)=Phi_elem(:,4)+Phi_e(:,i_onde)*exp(-1j*k_e(i_onde)*(n_phi'*(coord_elem(:,4)-x_centre)))*q(i_onde+ondes_element(ie)*(i_thetaphi-1));
                        
                    end
                end
                
                
                y=[coord_elem(2,:)];
                
                figure(10002)
                hold on
                c=[Phi_elem(8,1);Phi_elem(8,2);Phi_elem(8,3);Phi_elem(8,4)];
                plot(y,abs(c),'r+')
                
                figure(10001)
                hold on
                c=[Phi_elem(2,1);Phi_elem(2,2);Phi_elem(2,3);Phi_elem(2,4)]/(j*omega);
                plot(y,abs(c),'r+')
                
        end
        
    end
    
end



