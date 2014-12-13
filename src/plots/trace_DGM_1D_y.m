

for ie=1:nb_elements
    coord_elem=vcor(kconec(ie,:),1:2)';
    x_centre=centre_element(ie);

    e_edge=ie;
    parameter_element

    
    q=X_DGM_fortran(dof_start_element(ie)+[0:ondes_element(ie)*nb_theta-1]);
    Phi_elem=zeros(3,3);
    for i_thetaphi=1:nb_theta
        theta_phi=vec_theta(i_thetaphi);
        n_phi=[cos(theta_phi)*tau_x;sin(theta_phi)*tau_y];
        Phi_e=[cos(theta_phi);sin(theta_phi);-Z_e];
        Phi_elem(:,1)=Phi_elem(:,1)+Phi_e*exp(-j*k_e*(n_phi.'*(coord_elem(:,1)-x_centre)))*q(1+ondes_element(ie)*(i_thetaphi-1));
        Phi_elem(:,2)=Phi_elem(:,2)+Phi_e*exp(-j*k_e*(n_phi.'*(coord_elem(:,2)-x_centre)))*q(1+ondes_element(ie)*(i_thetaphi-1));
        Phi_elem(:,3)=Phi_elem(:,3)+Phi_e*exp(-j*k_e*(n_phi.'*(coord_elem(:,3)-x_centre)))*q(1+ondes_element(ie)*(i_thetaphi-1));
    end

    x=[coord_elem(1,:)];
    y=[coord_elem(2,:)];
    c=[Phi_elem(1,1);Phi_elem(1,2);Phi_elem(1,3)];
    c_trace=(c);

%     figure(1111)
%     hold on
%     patch(x,y,abs((c)));
%    colorbar


        % Trac? des figures pour le cas 1D
        figure(1000)
        for i_fig=2
            %subplot(2,2,i_fig)
            hold on
            plot(y,abs([Phi_elem(i_fig,1);Phi_elem(i_fig,2);Phi_elem(i_fig,3)]),'.')
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