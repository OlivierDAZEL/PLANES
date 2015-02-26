
coord_edge(1:2,1)=nodes(dirichlets(ie,1),1:2)';
coord_edge(1:2,2)=nodes(dirichlets(ie,2),1:2)';

e_edge=dirichlets(ie,3);
c_edge=mean(nodes(elements(e_edge,:),1:2))';


parameter_element

a=coord_edge(:,1);
b=coord_edge(:,2);


centre_temp=c_edge;
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

M_plus =B_e*M_e*W_e_plus *Lambda_e_plus;
M_moins=B_e*M_e*W_e_moins*Lambda_e_moins;


R_e=-inv(M_moins)*M_plus;


%F_e=A_x_fluid*nx+A_y_fluid*ny;
F_e=M_e*(Lambda_e_plus*W_e_plus+Lambda_e_moins*W_e_moins*R_e)*Omega_e_plus;


for i_thetapsi=1:nb_theta
    theta_psi=vec_theta(i_thetapsi);
    n_psi=[cos(theta_psi);sin(theta_psi)];
    Psi_e=conj(Phi_fluid(cos(theta_psi),sin(theta_psi),Z_e));
    ii=indice_fluid(e_edge,i_thetapsi,dof_start_element);

    for i_thetaphi=(1:nb_theta)
        theta_phi=vec_theta(i_thetaphi);
        n_phi=[cos(theta_phi);sin(theta_phi)];
        Phi_e=     Phi_fluid(cos(theta_phi),sin(theta_phi),Z_e) ;
        jj=indice_fluid(e_edge,i_thetaphi,dof_start_element);
        A(ii,jj)=A(ii,jj)+Psi_e'*F_e*Phi_e*...
            int_edge_2(j*k_e*n_psi,-j*k_e*n_phi,a,b,[c_edge c_edge]);
    end
end




% nx=cos(vec_theta);
% ny=sin(vec_theta);
% Phi=Phi_fluid_vector(nx,ny,Z_e,Shift_fluid);
% 
% 
% II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a,b,[c_edge c_edge]);
% MM=kron(II,-F_e);
% indice_test  =((1:nb_theta)-1)+dof_start_element(e_edge);  
% indice_champs=((1:nb_theta)-1)+dof_start_element(e_edge);
% A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;





