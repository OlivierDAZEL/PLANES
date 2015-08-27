% e'=e_1
% e=e_2
%

e_edge=e_1;
c_1=centre_element(e_1,nodes,elements);

parameter_element

c_e_p=c_e;
k_e_p=k_e;
Z_e_p=Z_e;
M_e_p=M_e;

e_edge=e_2;
c_2=centre_element(e_2,nodes,elements);
parameter_element


B_e=[nx ny 0; 0 0 1];
B_e_p=B_e;

W_e_plus= Phi_fluid( nx, ny,Z_e);
W_e_moins=Phi_fluid(-nx,-ny,Z_e);
W_e_0=Phi_fluid_0(nx,ny);


Omega_e=inv([W_e_plus W_e_moins W_e_0]);
Omega_e_plus=Omega_e(1,:);
Omega_e_moins=Omega_e(2,:);

Lambda_e_moins=omega/k_e;
Lambda_e_plus=-omega/k_e;

W_e_p_plus= Phi_fluid(-nx,-ny,Z_e_p);
W_e_p_moins=Phi_fluid( nx, ny,Z_e_p);
W_e_p_0=Phi_fluid_0(nx,ny);

Omega_e_p=inv([W_e_p_plus W_e_p_moins W_e_p_0]);
Omega_e_p_plus=Omega_e_p(1,:);
Omega_e_p_moins=Omega_e_p(2,:);


Lambda_e_p_moins=-omega/k_e_p;
Lambda_e_p_plus=  omega/k_e_p;

M1=[-B_e*M_e*W_e_plus*Lambda_e_plus B_e_p*M_e_p*W_e_p_plus*Lambda_e_p_plus];
M2=[ B_e*M_e*W_e_moins*Lambda_e_moins -B_e_p*M_e_p*W_e_p_moins*Lambda_e_p_moins];

Refl=inv(M2)*M1;

R_ee=Refl(1,1);
R_eep=Refl(1,2);
R_epe=Refl(2,1);
R_epep=Refl(2,2);


F_ee=  M_e*(W_e_plus*Lambda_e_plus+W_e_moins*Lambda_e_moins*R_ee)*Omega_e_plus;
F_eep= M_e*W_e_moins*Lambda_e_moins*R_eep*Omega_e_p_plus;
F_epe=-M_e_p*W_e_p_moins*Lambda_e_p_moins*R_epe*Omega_e_plus;
F_epep=-M_e_p*(W_e_p_moins*Lambda_e_p_moins*R_epep+W_e_p_plus*Lambda_e_p_plus)*Omega_e_p_plus;


for i_theta_test=1:nb_theta
    theta_test=vec_theta(i_theta_test);
    n_psi=[cos(theta_test);sin(theta_test)];
    for  i_theta_champs=(1:nb_theta)
        theta_champs=vec_theta(i_theta_champs);
        n_phi=[cos(theta_champs);sin(theta_champs)];
        % e'=e_1
        % e=e_2
        
        % Champs direct et test
        Psi_e=conj(Phi_fluid(cos(theta_test),sin(theta_test),Z_e)); 
        Phi_e=     Phi_fluid(cos(theta_champs),sin(theta_champs),Z_e);
        Psi_e_p=conj(Phi_fluid(cos(theta_test),sin(theta_test),Z_e_p)); 
        Phi_e_p=     Phi_fluid(cos(theta_champs),sin(theta_champs),Z_e_p);

        
        % ve^H F11 ue 
        indice_test  =indice_fluid(e_2,i_theta_test,dof_start_element);
        indice_champs=indice_fluid(e_2,i_theta_champs,dof_start_element);
        
        A(indice_test,indice_champs)=A(indice_test,indice_champs)+Psi_e'*F_ee *Phi_e...
            *int_edge_2(j*k_e*n_psi,-j*k_e*n_phi,a,b,[c_2 c_2]);
        
        % ve^H Feep ue'
        
        indice_test  =indice_fluid(e_2,i_theta_test,dof_start_element);
        indice_champs=indice_fluid(e_1,i_theta_champs,dof_start_element);

        A(indice_test,indice_champs)=A(indice_test,indice_champs)+Psi_e'*F_eep *Phi_e_p...
            *int_edge_2(j*k_e*n_psi,-j*k_e_p*n_phi,a,b,[c_2 c_1]);
        % -ve'^H F+ ue
        
        indice_test  =indice_fluid(e_1,i_theta_test,dof_start_element);
        indice_champs=indice_fluid(e_2,i_theta_champs,dof_start_element);
        
       % F_epe
        
        
        A(indice_test,indice_champs)=A(indice_test,indice_champs)+Psi_e_p'*F_epe *Phi_e...
            *int_edge_2(j*k_e_p*n_psi,-j*k_e*n_phi,a,b,[c_1 c_2]);
        % -ve'^H F- ue'
        
        indice_test  =indice_fluid(e_1,i_theta_test,dof_start_element);
        indice_champs=indice_fluid(e_1,i_theta_champs,dof_start_element);
        A(indice_test,indice_champs)=A(indice_test,indice_champs)+Psi_e_p'*F_epep *Phi_e_p...
            *int_edge_2(j*k_e_p*n_psi,-j*k_e_p*n_phi,a,b,[c_1 c_1]);
    end
end