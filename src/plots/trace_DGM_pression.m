% trace_DGM_pression.m
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

    vertices=[vcor(kconec(ie,:),1:2)';0 0 0];
    faces = [1 2 3]';

    for ii=1:2
        [vertices, faces]=linearSubdivision(vertices, faces);
    end
    x_centre=centre_element(ie);

    e_edge=ie;
    parameter_element
    x_centre=diag([tau_x tau_y])*x_centre;
    
    q=X_DGM(dof_start_element(ie)+[0:ondes_element(ie)*nb_theta-1]);

    vec_normales=[cos(vec_theta);sin(vec_theta)]';

    if ((element_label(ie)==0)|(floor(element_label(ie)/1000)==2)|(floor(element_label(ie)/1000)==3)|(floor(element_label(ie)/1000)==8))
        for i_point=1:size(vertices,2);
            x_point=diag([tau_x tau_y])*vertices(1:2,i_point);
            Delta=diag(exp(-j*k_e*vec_normales*(x_point-x_centre)));
            Phi_point(1:3,i_point)=sum(Phi_fluid(cos(vec_theta),sin(vec_theta),Z_e)*Delta*q,2);
        end
        figure(1111)
        for i_faces=1:size(faces,2)
            patch(vertices(1,faces(:,i_faces)),vertices(2,faces(:,i_faces)),abs(Phi_point(3,faces(:,i_faces))'));
        end
    elseif (floor(element_label(e_edge)/1000)==4)
        for i_point=1:size(vertices,2);
            x_point=vertices(1:2,i_point);
            Phi_plot=Phi_Biot(cos(vec_theta),sin(vec_theta),delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);

            Delta=diag(exp(-j*delta_1*vec_normales*(x_point-x_centre)));
            Phi_point(1:8,i_point)=sum(Phi_plot(1:8,1:nb_theta)*Delta*q(1:3:end),2);

            Delta=diag(exp(-j*delta_2*vec_normales*(x_point-x_centre)));
            Phi_point(1:8,i_point)=Phi_point(1:8,i_point)+sum(Phi_plot(1:8,nb_theta+(1:nb_theta))*Delta*q(2:3:end),2);

            Delta=diag(exp(-j*delta_3*vec_normales*(x_point-x_centre)));
            Phi_point(1:8,i_point)=Phi_point(1:8,i_point)+sum(Phi_plot(1:8,2*nb_theta+(1:nb_theta))*Delta*q(3:3:end),2);
        end


        figure(1111)

        for i_faces=1:size(faces,2)
            patch(vertices(1,faces(:,i_faces)),vertices(2,faces(:,i_faces)),abs(Phi_point(8,faces(:,i_faces))'));
        end

        figure(2222)
        for i_faces=1:size(faces,2)
            patch(vertices(1,faces(:,i_faces)),vertices(2,faces(:,i_faces)),sqrt((Phi_point(1,faces(:,i_faces))').^2+(Phi_point(2,faces(:,i_faces))').^2));
        end        
        


    end



end

axis equal
shading interp

figure(1111)
colorbar

for ie=1:nb_elements
    vertices=[vcor(kconec(ie,:),1:2)];
    faces = [1 2 3]';

    line([vertices(faces(1),1) vertices(faces(2),1)],[vertices(faces(1),2) vertices(faces(2),2)],'Linewidth',1,'Color','k')
    line([vertices(faces(2),1) vertices(faces(3),1)],[vertices(faces(2),2) vertices(faces(3),2)],'Linewidth',1,'Color','k')
    line([vertices(faces(3),1) vertices(faces(1),1)],[vertices(faces(3),2) vertices(faces(1),2)],'Linewidth',1,'Color','k')

end



% if savefigure
%     imprime_figure_DGM
% end