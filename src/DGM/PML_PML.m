% e'=e_1
% e=e_2
%

e_edge=e_1;
c_1=centre_element(e_1,nodes,elements);

parameter_element

tau_x_e1=tau_x;
tau_y_e1=tau_y;
Z_e1=Z_e;

c_e_p=c_e;
k_e_p=k_e;
Z_e_p=Z_e;
M_e_p=M_e;

e_edge=e_2;
c_2=centre_element(e_2,nodes,elements);
parameter_element

Z_e2=Z_e;

tau_x_e2=tau_x;
tau_y_e2=tau_y;

	
ne1=sqrt((nx/tau_x_e1)^2+(ny/tau_y_e1)^2);
ne2=sqrt((nx/tau_x_e2)^2+(ny/tau_y_e2)^2);


unc=1;
deuxc=2;


Re1e1=(Z_e2*tau_x_e2*tau_y_e2*ne2-Z_e1*tau_x_e1*tau_y_e1*ne1)/(Z_e2*tau_x_e2*tau_y_e2*ne2+Z_e1*tau_x_e1*tau_y_e1*ne1);
Re1e2=(2.00D+00)*Z_e2*tau_x_e2*tau_y_e2*ne2^2/(Z_e1*ne1*tau_x_e2*tau_y_e2*ne2+Z_e2*tau_x_e1*tau_y_e1*ne1^2);
Re2e1=(2.00D+00)*Z_e1*tau_x_e1*tau_y_e1*ne1^2/(Z_e1*tau_x_e2*tau_y_e2*ne2^2+Z_e2*ne2*tau_x_e1*tau_y_e1*ne1);
Re2e2=(Z_e1*tau_x_e1*tau_y_e1*ne1-Z_e2*tau_x_e2*tau_y_e2*ne2)/(Z_e2*tau_x_e2*tau_y_e2*ne2+Z_e1*tau_x_e1*tau_y_e1*ne1);
 
 
F_e1e1(1,1)=-Z_e1*(nx^2)*(unc+Re1e1)/(deuxc*tau_x_e1^2*ne1);
F_e1e1(1,2)=-Z_e1*(nx*ny)*(unc+Re1e1)/(deuxc*tau_x_e1*tau_y_e1*ne1);
F_e1e1(1,3)=-(nx)*(unc+Re1e1)/(deuxc*tau_x_e1);
F_e1e1(2,1)=-Z_e1*(nx*ny)*(unc+Re1e1)/(deuxc*tau_x_e1*tau_y_e1*ne1);
F_e1e1(2,2)=-Z_e1*(ny^2)*(unc+Re1e1)/(deuxc*tau_y_e1^2*ne1);
F_e1e1(2,3)=-(ny)*(unc+Re1e1)/(deuxc*tau_y_e1);
F_e1e1(3,1)=(nx)*(Re1e1-unc)/(deuxc*tau_x_e1);
F_e1e1(3,2)=(ny)*(Re1e1-unc)/(deuxc*tau_y_e1);
F_e1e1(3,3)=ne1*(Re1e1-unc)/(deuxc*Z_e1);
 
 
F_e2e2(1,1)=-Z_e2*(nx^2)*(unc+Re2e2)/(deuxc*tau_x_e2^2*ne2);
F_e2e2(1,2)=-Z_e2*(nx*ny)*(unc+Re2e2)/(deuxc*tau_x_e2*tau_y_e2*ne2);
F_e2e2(1,3)=(nx)*(unc+Re2e2)/(deuxc*tau_x_e2);
F_e2e2(2,1)=-Z_e2*(nx*ny)*(unc+Re2e2)/(deuxc*tau_x_e2*tau_y_e2*ne2);
F_e2e2(2,2)=-Z_e2*(ny^2)*(unc+Re2e2)/(deuxc*tau_y_e2^2*ne2);
F_e2e2(2,3)=(ny)*(unc+Re2e2)/(deuxc*tau_y_e2);
F_e2e2(3,1)=-(nx)*(Re2e2-unc)/(deuxc*tau_x_e2);
F_e2e2(3,2)=-(ny)*(Re2e2-unc)/(deuxc*tau_y_e2);
F_e2e2(3,3)=ne2*(Re2e2-unc)/(deuxc*Z_e2);
 
 
F_e2e1(1,1)=Z_e2*(nx^2)/(deuxc*tau_x_e2*tau_x_e1*ne1^2);
F_e2e1(1,2)=Z_e2*(nx*ny)/(deuxc*tau_x_e2*tau_y_e1*ne1^2);
F_e2e1(1,3)=(nx)*Z_e2/(deuxc*Z_e1*tau_x_e2*ne1);
F_e2e1(2,1)=Z_e2*(nx*ny)/(deuxc*tau_x_e1*tau_y_e2*ne1^2);
F_e2e1(2,2)=Z_e2*(ny^2)/(deuxc*tau_y_e2*tau_y_e1*ne1^2);
F_e2e1(2,3)=(ny)*Z_e2/(deuxc*Z_e1*tau_y_e2*ne1);
F_e2e1(3,1)=ne2*(nx)/(deuxc*tau_x_e1*ne1^2);
F_e2e1(3,2)=ne2*(ny)/(deuxc*tau_y_e1*ne1^2);
F_e2e1(3,3)=ne2*unc/(deuxc*Z_e1*ne1);
 
F_e2e1=Re2e1*F_e2e1*ne2;
 
 
F_e1e2(1,1)=Z_e1*(nx^2)/(deuxc*tau_x_e1*tau_x_e2*ne2^2);
F_e1e2(1,2)=Z_e1*(nx*ny)/(deuxc*tau_x_e1*tau_y_e2*ne2^2);
F_e1e2(1,3)=-(nx)*Z_e1/(deuxc*Z_e2*tau_x_e1*ne2);
F_e1e2(2,1)=Z_e1*(nx*ny)/(deuxc*tau_x_e2*tau_y_e1*ne2^2);
F_e1e2(2,2)=Z_e1*(ny^2)/(deuxc*tau_y_e1*tau_y_e2*ne2^2);
F_e1e2(2,3)=-(ny)*Z_e1/(deuxc*Z_e2*tau_y_e1*ne2);
F_e1e2(3,1)=-ne1*(nx)/(deuxc*tau_x_e2*ne2^2);
F_e1e2(3,2)=-ne1*(ny)/(deuxc*tau_y_e2*ne2^2);
F_e1e2(3,3)=ne1*unc/(deuxc*Z_e2*ne2);
 
F_e1e2=Re1e2*F_e1e2*ne1;



nx=cos(vec_theta);
ny=sin(vec_theta);
Phi=Phi_fluid_vector(nx,ny,air.Z,Shift_fluid);

II=int_edge_2vectorielle(1j*delta_test*[nx*tau_x_e2;ny*tau_y_e2],-1j*delta_champs*[nx*tau_x_e2;ny*tau_y_e2],a,b,[c_2 c_2]);
MM=kron(II,F_e2e2);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_2);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_2);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi*tau_x_e2*tau_y_e2;

        
II=int_edge_2vectorielle(1j*delta_test*[nx*tau_x_e2;ny*tau_y_e2],-1j*delta_champs*[nx*tau_x_e1;ny*tau_y_e1],a,b,[c_2 c_1]);
MM=kron(II,F_e2e1);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_2);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_1);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi*tau_x_e2*tau_y_e2;

II=int_edge_2vectorielle(1j*delta_test*[nx*tau_x_e1;ny*tau_y_e1],-1j*delta_champs*[nx*tau_x_e2;ny*tau_y_e2],a,b,[c_1 c_2]);
MM=kron(II,F_e1e2);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_1);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_2);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi*tau_x_e1*tau_y_e1;

II=int_edge_2vectorielle(1j*delta_test*[nx*tau_x_e1;ny*tau_y_e1],-1j*delta_champs*[nx*tau_x_e1;ny*tau_y_e1],a,b,[c_1 c_1]);
MM=kron(II,F_e1e1);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_1);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_1);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi*tau_x_e1*tau_y_e1;

