e_edge=e_1;
parameter_element


[F_plus,F_moins]=Split_fluid(nx,ny,k_e,Z_e,omega,M_e);

delta_test=k_e;
delta_champs=k_e;

nx=cos(vec_theta);
ny=sin(vec_theta);
Phi=Phi_fluid_vector(nx,ny,Z_e,Shift_fluid);

II=int_edge_2vectorielle(1j*delta_test*[nx;ny],-1j*delta_champs*[nx;ny],a,b,[c_2 c_2]);
MM=kron(II,F_plus);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_2);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_2);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;

        
II=int_edge_2vectorielle(1j*delta_test*[nx;ny],-1j*delta_champs*[nx;ny],a,b,[c_2 c_1]);
MM=kron(II,F_moins);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_2);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_1);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;

II=int_edge_2vectorielle(1j*delta_test*[nx;ny],-1j*delta_champs*[nx;ny],a,b,[c_1 c_2]);
MM=kron(II,-F_plus);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_1);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_2);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;

II=int_edge_2vectorielle(1j*delta_test*[nx;ny],-1j*delta_champs*[nx;ny],a,b,[c_1 c_1]);
MM=kron(II,-F_moins);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_1);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_1);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;

