parameter_element

W_e_plus= Phi_Biot(-nx,-ny,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
W_e_moins=Phi_Biot( nx, ny,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
W_e_0=Phi_Biot_0(nx,ny);

Omega_e=inv([W_e_plus W_e_moins W_e_0]);

Omega_e_plus=Omega_e(1:3,:);
Omega_e_moins=Omega_e(4:6,:);

Lambda_e_moins=-diag([omega/delta_1 omega/delta_2 omega/delta_3]);
Lambda_e_plus=-Lambda_e_moins;

F_plus=M_e*W_e_plus*Lambda_e_plus*Omega_e_plus;
F_moins=M_e*W_e_moins*Lambda_e_moins*Omega_e_moins;

%F_plus+F_moins



delta=[delta_1 delta_2 delta_3];


for i_thetapsi=1:nb_theta
    theta_test=vec_theta(i_thetapsi);
    n_psi=[cos(theta_test);sin(theta_test)];
    Psi_e=conj(Phi_Biot(cos(theta_test)  ,sin(theta_test)  ,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega));
    for  i_thetaphi=(1:nb_theta)
        theta_champs=vec_theta(i_thetaphi);
        n_phi=[cos(theta_champs);sin(theta_champs)];

        Phi_e=     Phi_Biot(cos(theta_champs),sin(theta_champs),delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
        i_test=1:3; % Balayage des ondes de Biot champs test
        i_champs=1:3;

        
        integrale_border11=int_edge_2vectorielle(j*(ones(2,1)*delta).*(n_psi*ones(1,3)),-j*(ones(2,1)*delta).*(n_phi*ones(1,3)),a,b,[c_1 c_1]);
        integrale_border12=int_edge_2vectorielle(j*(ones(2,1)*delta).*(n_psi*ones(1,3)),-j*(ones(2,1)*delta).*(n_phi*ones(1,3)),a,b,[c_1 c_2]);
        integrale_border21=int_edge_2vectorielle(j*(ones(2,1)*delta).*(n_psi*ones(1,3)),-j*(ones(2,1)*delta).*(n_phi*ones(1,3)),a,b,[c_2 c_1]);
        integrale_border22=int_edge_2vectorielle(j*(ones(2,1)*delta).*(n_psi*ones(1,3)),-j*(ones(2,1)*delta).*(n_phi*ones(1,3)),a,b,[c_2 c_2]);

        
        %ve^H F+ ue
        ii=indice_Biot(e_2,i_thetapsi,i_test,dof_start_element);
        jj=indice_Biot(e_2,i_thetaphi,i_champs,dof_start_element);
        A(ii,jj)=A(ii,jj)+(Psi_e(:,i_test)'*F_plus *Phi_e(:,i_champs)).*integrale_border22;
        %ve^H F- ue'
        ii=indice_Biot(e_2,i_thetapsi,i_test,dof_start_element);
        jj=indice_Biot(e_1,i_thetaphi,i_champs,dof_start_element);
        A(ii,jj)=A(ii,jj)+Psi_e(:,i_test)'*F_moins *Phi_e(:,i_champs).*integrale_border21;
        %-ve'^H F+ ue
        ii=indice_Biot(e_1,i_thetapsi,i_test,dof_start_element);
        jj=indice_Biot(e_2,i_thetaphi,i_champs,dof_start_element);
        A(ii,jj)=A(ii,jj)-Psi_e(:,i_test)'*F_plus *Phi_e(:,i_champs).*integrale_border12;
        %-ve'^H F- ue'
        ii=indice_Biot(e_1,i_thetapsi,i_test,dof_start_element);
        jj=indice_Biot(e_1,i_thetaphi,i_champs,dof_start_element);
        A(ii,jj)=A(ii,jj)-Psi_e(:,i_test)'*F_moins *Phi_e(:,i_champs).*integrale_border11;




    end
end




% for i_thetapsi=1:nb_theta
%     theta_test=vec_theta(i_thetapsi);
%     n_psi=[cos(theta_test);sin(theta_test)];
%     for  i_thetaphi=(1:nb_theta)
%         theta_champs=vec_theta(i_thetaphi);
%         n_phi=[cos(theta_champs);sin(theta_champs)];
%         Psi_e=conj(Phi_Biot(cos(theta_test)  ,sin(theta_test)  ,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega));
%         Phi_e=     Phi_Biot(cos(theta_champs),sin(theta_champs),delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega);
%         
%         integrale_border22=int_edge_2vectorielle(j*(ones(2,1)*delta).*(n_psi*ones(1,3)),-j*(ones(2,1)*delta).*(n_phi*ones(1,3)),a,b,[c_2 c_2])
%         
%         
%         for i_test=1:3 % Balayage des ondes de Biot champs test
%             for i_champs=1:3 % Balayage des ondes de Biot champs inconnu
%                 % ve^H F+ ue
%                 ii=indice_Biot(e_2,i_thetapsi,i_test);
%                 jj=indice_Biot(e_2,i_thetaphi,i_champs);
%                 M_DGM(ii,jj)=M_DGM(ii,jj)+Psi_e(:,i_test)'*F_plus *Phi_e(:,i_champs)...
%                         *int_edge_2(j*delta(i_test)*n_psi,-j*delta(i_champs)*n_phi,a,b,[c_2 c_2]);
%                 % ve^H F- ue'
%                 ii=indice_Biot(e_2,i_thetapsi,i_test);
%                 jj=indice_Biot(e_1,i_thetaphi,i_champs);
%                 M_DGM(ii,jj)=M_DGM(ii,jj)+Psi_e(:,i_test)'*F_moins *Phi_e(:,i_champs)...
%                         *int_edge_2(j*delta(i_test)*n_psi,-j*delta(i_champs)*n_phi,a,b,[c_2 c_1]);
%                 % -ve'^H F+ ue
%                 ii=indice_Biot(e_1,i_thetapsi,i_test);
%                 jj=indice_Biot(e_2,i_thetaphi,i_champs);
%                 M_DGM(ii,jj)=M_DGM(ii,jj)-Psi_e(:,i_test)'*F_plus *Phi_e(:,i_champs)...
%                         *int_edge_2(j*delta(i_test)*n_psi,-j*delta(i_champs)*n_phi,a,b,[c_1 c_2]);
%                 % -ve'^H F- ue'
%                 ii=indice_Biot(e_1,i_thetapsi,i_test);
%                 jj=indice_Biot(e_1,i_thetaphi,i_champs);
%                 M_DGM(ii,jj)=M_DGM(ii,jj)-Psi_e(:,i_test)'*F_moins *Phi_e(:,i_champs)...
%                     *int_edge_2(j*delta(i_test)*n_psi,-j*delta(i_champs)*n_phi,a,b,[c_1 c_1]);
% 
%             end
%         end
%     end
% end