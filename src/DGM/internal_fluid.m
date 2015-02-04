c_e=air.c;
k_e=omega/c_e;
Z_e=air.Z;

M_e=diag([air.rho,air.rho,1/air.K]);

[F_plus,F_moins]=Split_fluid(nx,ny,k_e,Z_e,omega,M_e);

delta_test=k_e;
delta_champs=k_e;


for i_theta_test=1:nb_theta
    theta_test=vec_theta(i_theta_test);
    n_psi=[cos(theta_test);sin(theta_test)];
    for  i_theta_champs=(1:nb_theta)
        theta_champs=vec_theta(i_theta_champs);
        n_phi=[cos(theta_champs);sin(theta_champs)];
        
        % Champs direct et test
        Psi_e=conj(Phi_fluid(cos(theta_test),sin(theta_test),Z_e));
        Phi_e=     Phi_fluid(cos(theta_champs),sin(theta_champs),Z_e);
        
        indice_test  =(i_theta_test-1)+dof_start_element(e_2);  
        indice_champs=(i_theta_champs-1)+dof_start_element(e_2);
        
        A(indice_test,indice_champs)=A(indice_test,indice_champs)+Psi_e'*F_plus *Phi_e...
            *int_edge_2(j*delta_test*n_psi,-j*delta_champs*n_phi,a,b,[c_2 c_2]);
        % ve^H F- ue'
        indice_test  =(i_theta_test-1)+dof_start_element(e_2);  
        indice_champs=(i_theta_champs-1)+dof_start_element(e_1);
        
        A(indice_test,indice_champs)=A(indice_test,indice_champs)+Psi_e'*F_moins *Phi_e...
            *int_edge_2(j*delta_test*n_psi,-j*delta_champs*n_phi,a,b,[c_2 c_1]);
        % -ve'^H F+ ue
        indice_test  =(i_theta_test-1)  +dof_start_element(e_1);  
        indice_champs=(i_theta_champs-1)+dof_start_element(e_2);
        A(indice_test,indice_champs)=A(indice_test,indice_champs)-Psi_e'*F_plus *Phi_e...
            *int_edge_2(j*delta_test*n_psi,-j*delta_champs*n_phi,a,b,[c_1 c_2]);
        % -ve'^H F- ue'
        indice_test  =(i_theta_test-1)  +dof_start_element(e_1);  
        indice_champs=(i_theta_champs-1)+dof_start_element(e_1);
        A(indice_test,indice_champs)=A(indice_test,indice_champs)-Psi_e'*F_moins *Phi_e...
            *int_edge_2(j*delta_test*n_psi,-j*delta_champs*n_phi,a,b,[c_1 c_1]);
    end
end

