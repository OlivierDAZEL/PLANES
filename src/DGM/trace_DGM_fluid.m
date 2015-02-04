
figure
hold on



for ie=1:nb_elements
    coord_elem=nodes(elements(ie,:),1:2)';
    x_centre=mean(nodes(elements(ie,:),1:2))';

    e_edge=ie;
    parameter_element
    
    q=X(dof_start_element(ie)+[0:ondes_element(ie)*nb_theta-1]);
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
    c=[Phi_elem(3,1);Phi_elem(3,2);Phi_elem(3,3)];
    c_trace=(c);

    plot(x,abs(c),'r+')

end