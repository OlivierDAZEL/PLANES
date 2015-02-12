% e'=e_1
% e=e_2


e_edge=e_1;
parameter_element

M_e_p=M_e;

W_e_p_plus= Phi_fluid( nx, ny,Z_e);
W_e_p_moins=Phi_fluid(-nx,-ny,Z_e);
W_e_p_0=Phi_fluid_0(nx,ny);


Omega_e_p=inv([W_e_p_plus W_e_p_moins W_e_p_0]);
Omega_e_p_plus=Omega_e_p(1,:);
Omega_e_p_moins=Omega_e_p(2,:);

Lambda_e_p_moins=omega/k_e;
Lambda_e_p_plus=-Lambda_e_p_moins;

k_e_p=k_e;
Z_e_p=Z_e;

e_edge=e_2;
parameter_element


W_e_plus= Phi_Biot(-nx,-ny,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
W_e_moins=Phi_Biot( nx, ny,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
W_e_0=Phi_Biot_0(nx,ny);

Omega_e=inv([W_e_plus W_e_moins W_e_0]);
Omega_e_plus=Omega_e(1:3,:);
Omega_e_moins=Omega_e(4:6,:);

Lambda_e_moins=-diag([omega/delta_1 omega/delta_2 omega/delta_3]);
Lambda_e_plus=-Lambda_e_moins;

k_e=[delta_1 delta_2 delta_3];


B_e_p=zeros(4,3);
B_e=zeros(4,8);

B_e_p(3,1)=nx;
B_e_p(3,2)=ny;
B_e_p(4,3)=1;

B_e(1,1)=1;
B_e(2,2)=1;
B_e(3,3)=nx;
B_e(3,4)=ny;
B_e(4,8)=1;


M1=[-B_e*M_e*W_e_plus*Lambda_e_plus B_e_p*M_e_p*W_e_p_plus*Lambda_e_p_plus];
M2=[ B_e*M_e*W_e_moins*Lambda_e_moins -B_e_p*M_e_p*W_e_p_moins*Lambda_e_p_moins];

Refl=inv(M2)*M1;

R_ee=Refl(1:3,1:3);
R_eep=Refl(1:3,4);
R_epe=Refl(4,1:3);
R_epep=Refl(4,4);


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
        Psi_e_p=conj(Phi_fluid(cos(theta_test),sin(theta_test),Z_e_p));
        Phi_e_p=     Phi_fluid(cos(theta_champs),sin(theta_champs),Z_e_p);
        Psi_e=conj(Phi_Biot(cos(theta_test),sin(theta_test),delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega));
        Phi_e=     Phi_Biot(cos(theta_champs),sin(theta_champs),delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
     % ve^H F11 ue

        for i_test=1:3
            for i_champs=1:3
                indice_test  =indice_Biot(e_2,i_theta_test,i_test,dof_start_element);
                indice_champs=indice_Biot(e_2,i_theta_champs,i_champs,dof_start_element);

                A(indice_test,indice_champs)=A(indice_test,indice_champs)+Psi_e(:,i_test)'*F_ee *Phi_e(:,i_champs)...
                    *int_edge_2(j*k_e(i_test)*n_psi,-j*k_e(i_champs)*n_phi,a,b,[c_2 c_2]);

            end
        end
        % ve^H Feep ue'
        for i_test=1:3
            indice_test  =indice_Biot(e_2,i_theta_test,i_test,dof_start_element);
            indice_champs=indice_fluid(e_1,i_theta_champs,dof_start_element);

            A(indice_test,indice_champs)=A(indice_test,indice_champs)+Psi_e(:,i_test)'*F_eep *Phi_e_p...
                *int_edge_2(j*k_e(i_test)*n_psi,-j*k_e_p*n_phi,a,b,[c_2 c_1]);
        end

        % -ve'^H F+ ue
        for i_champs=1:3
            indice_test  =indice_fluid(e_1,i_theta_test,dof_start_element);
            indice_champs=indice_Biot(e_2,i_theta_champs,i_champs,dof_start_element);
            A(indice_test,indice_champs)=A(indice_test,indice_champs)+Psi_e_p'*F_epe *Phi_e(:,i_champs)...
                *int_edge_2(j*k_e_p*n_psi,-j*k_e(i_champs)*n_phi,a,b,[c_1 c_2]);
            % -ve'^H F- ue'
        end

        indice_test  =indice_fluid(e_1,i_theta_test,dof_start_element);
        indice_champs=indice_fluid(e_1,i_theta_champs,dof_start_element);
        A(indice_test,indice_champs)=A(indice_test,indice_champs)+Psi_e_p'*F_epep *Phi_e_p...
            *int_edge_2(j*k_e_p*n_psi,-j*k_e_p*n_phi,a,b,[c_1 c_1]);
    end
end