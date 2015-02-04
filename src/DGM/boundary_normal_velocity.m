
coord_edge(1:2,1)=nodes(loads(ie,1),1:2)';
coord_edge(1:2,2)=nodes(loads(ie,2),1:2)';


e_edge=loads(ie,3);
c_edge=mean(nodes(elements(e_edge,:),1:2))';
parameter_element


a=coord_edge(:,1);
b=coord_edge(:,2);

centre_temp=mean(nodes(elements(e_edge,:),1:2))';
centre_edge=(a+b)/2;

n_centre=centre_temp-centre_edge;

ne=normal_edge(coord_edge);

if (n_centre'*ne<0)
    ne=-ne;
end


nx=ne(1);
ny=ne(2);


W_e_plus=Phi_fluid(nx,ny,Z_e);
W_e_moins=Phi_fluid(-nx,-ny,Z_e);
W_e_0=Phi_fluid_0(nx,ny);


W_e=[W_e_plus W_e_moins W_e_0];
Omega_e=inv(W_e);

Omega_e_plus=Omega_e(1,:);
Omega_e_moins=Omega_e(2,:);


Lambda_e_plus= diag(c_e);
Lambda_e_moins=-Lambda_e_plus;

B_e=[0 0 1];
s_e=1;

M_plus =B_e*M_e*W_e_plus *Lambda_e_plus;
M_moins=B_e*M_e*W_e_moins*Lambda_e_moins;


R_e=-inv(M_moins)*M_plus;
S_e= inv(M_moins)*s_e;


%F_e=A_x_fluid*nx+A_y_fluid*ny;
F_e=M_e*(Lambda_e_plus*W_e_plus+Lambda_e_moins*W_e_moins*R_e)*Omega_e_plus;
S_e=M_e*                        Lambda_e_moins*W_e_moins*S_e;


for i_theta_test=1:nb_theta
    theta_psi=vec_theta(i_theta_test);
    n_psi=[cos(theta_psi);sin(theta_psi)];
    for i_theta_champs=(1:nb_theta)
        theta_phi=vec_theta(i_theta_champs);
        n_phi=[cos(theta_phi);sin(theta_phi)];
        Psi_e=conj(Phi_fluid(cos(theta_psi),sin(theta_psi),Z_e));
        Phi_e=     Phi_fluid(cos(theta_phi),sin(theta_phi),Z_e) ;
        
        indice_test  =(i_theta_test-1)  +dof_start_element(e_edge);  
        indice_champs=(i_theta_champs-1)+dof_start_element(e_edge);
        
        A(indice_test,indice_champs)=A(i_theta_test,indice_champs)+Psi_e'*F_e*Phi_e*...
            int_edge_2(j*k_e*n_psi,-j*k_e*n_phi,a,b,[c_edge c_edge]);
    end
end


for i_theta_test=1:nb_theta
    theta_psi=vec_theta(i_theta_test);
    n_psi=[cos(theta_psi);sin(theta_psi)];
    Psi_e=conj(Phi_fluid(cos(theta_psi),sin(theta_psi),Z_e));
    indice_test =(i_theta_test-1)+dof_start_element(e_edge);  
    F(indice_test)=F(indice_test)+Psi_e'*S_e*...
        int_edge_1(j*k_e*n_psi,a,b,c_edge);
end







