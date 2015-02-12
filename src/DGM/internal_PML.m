[tau_x,tau_y]=parameter_PML(element_label(e_edge));

[F_plus,F_moins]=Split_PML(nx,ny,Z_e,tau_x,tau_y);


nx=cos(vec_theta);
ny=sin(vec_theta);
Phi=Phi_fluid_vector(nx,ny,air.Z,Shift_fluid);

II=int_edge_2vectorielle(1j*k_air*[nx*tau_x;ny*tau_y],-1j*k_air*[nx*tau_x;ny*tau_y],a,b,[c_2 c_2]);
MM=kron(II,F_plus);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_2);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_2);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi*tau_x*tau_y;

        
II=int_edge_2vectorielle(1j*k_air*[nx*tau_x;ny*tau_y],-1j*k_air*[nx*tau_x;ny*tau_y],a,b,[c_2 c_1]);
MM=kron(II,F_moins);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_2);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_1);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi*tau_x*tau_y;

II=int_edge_2vectorielle(1j*k_air*[nx*tau_x;ny*tau_y],-1j*k_air*[nx*tau_x;ny*tau_y],a,b,[c_1 c_2]);
MM=kron(II,-F_plus);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_1);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_2);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi*tau_x*tau_y;

II=int_edge_2vectorielle(1j*k_air*[nx*tau_x;ny*tau_y],-1j*k_air*[nx*tau_x;ny*tau_y],a,b,[c_1 c_1]);
MM=kron(II,-F_moins);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_1);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_1);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi*tau_x*tau_y;

