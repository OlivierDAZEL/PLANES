


figure(10002)
hold on




for ie=1:nb_elements
    coord_elem=nodes(elements(ie,:),1:2)';
    x_centre=mean(nodes(elements(ie,:),1:2))';
    
    e_edge=ie;
    parameter_element
    
    q=X(dof_start_element(ie)+[0:ondes_element(ie)*nb_theta-1]);
    
    switch floor(element_label(ie)/1000)
        case {0,8}

            Phi_elem=zeros(3,3);
            for i_thetaphi=1:nb_theta
                theta_phi=vec_theta(i_thetaphi);
                n_phi=[cos(theta_phi)*tau_x;sin(theta_phi)*tau_y];
                Phi_e=Phi_fluid(cos(theta_phi),sin(theta_phi),Z_e);
                Phi_elem(:,1)=Phi_elem(:,1)+Phi_e*exp(-j*k_e*(n_phi.'*(coord_elem(:,1)-x_centre)))*q(1+ondes_element(ie)*(i_thetaphi-1));
                Phi_elem(:,2)=Phi_elem(:,2)+Phi_e*exp(-j*k_e*(n_phi.'*(coord_elem(:,2)-x_centre)))*q(1+ondes_element(ie)*(i_thetaphi-1));
                Phi_elem(:,3)=Phi_elem(:,3)+Phi_e*exp(-j*k_e*(n_phi.'*(coord_elem(:,3)-x_centre)))*q(1+ondes_element(ie)*(i_thetaphi-1));
            end
            
            x=[coord_elem(1,:)];
            y=[coord_elem(2,:)];

            figure(10002)
            hold on
            c=[Phi_elem(3,1);Phi_elem(3,2);Phi_elem(3,3)];
            plot(y,abs(c),'r+')
            
        case {4,5}    
            
k_e=[delta_1 delta_2 delta_3];
Phi_elem=zeros(8,3);
for i_thetaphi=1:nb_theta
    theta_phi=vec_theta(i_thetaphi);
    n_phi=[cos(theta_phi);sin(theta_phi)];
    Phi_e=Phi_Biot(cos(theta_phi),sin(theta_phi),delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
    for i_onde=1:3
        Phi_elem(:,1)=Phi_elem(:,1)+Phi_e(:,i_onde)*exp(-j*k_e(i_onde)*(n_phi'*(coord_elem(:,1)-x_centre)))*q(i_onde+ondes_element(ie)*(i_thetaphi-1));
        Phi_elem(:,2)=Phi_elem(:,2)+Phi_e(:,i_onde)*exp(-j*k_e(i_onde)*(n_phi'*(coord_elem(:,2)-x_centre)))*q(i_onde+ondes_element(ie)*(i_thetaphi-1));
        Phi_elem(:,3)=Phi_elem(:,3)+Phi_e(:,i_onde)*exp(-j*k_e(i_onde)*(n_phi'*(coord_elem(:,3)-x_centre)))*q(i_onde+ondes_element(ie)*(i_thetaphi-1));
    end
end

            
            x=[coord_elem(1,:)];
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