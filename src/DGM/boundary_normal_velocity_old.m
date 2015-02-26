%%%%% Coordinates of the edge's vertices
coord_edge(1:2,1)=nodes(loads(ie,1),1:2)';
coord_edge(1:2,2)=nodes(loads(ie,2),1:2)';

a=coord_edge(:,1);
b=coord_edge(:,2);

h=norm(b-a);
n_ab=(b-a)/h;

%%%%% Element linked to the edge

e_edge=loads(ie,3);
c_edge=mean(nodes(elements(e_edge,:),1:2))';

parameter_element

%%%%% vector normal pointing towards \Omega_e

centre_temp=mean(nodes(elements(e_edge,:),1:2))';
centre_edge=(a+b)/2;
n_centre=centre_temp-centre_edge;
ne=normal_edge(coord_edge);
if (n_centre'*ne<0)
    ne=-ne;
end
nx=ne(1);
ny=ne(2);



F_e=zeros(3,3);
S_e=zeros(3,1);
F_e(1,1)=-Z_e*(nx^2);
F_e(1,2)=-Z_e*(nx*ny);
F_e(1,3)= nx;
F_e(2,1)=-Z_e*(nx*ny);
F_e(2,2)=-Z_e*(ny^2);
F_e(2,3)= ny;

S_e(1)=1j*omega*Z_e*nx;
S_e(2)=1j*omega*Z_e*ny;
S_e(3)=1j*omega;



% 
% 
% W_e_plus=Phi_fluid(nx,ny,Z_e);
% W_e_moins=Phi_fluid(-nx,-ny,Z_e);
% W_e_0=Phi_fluid_0(nx,ny);
% 
% 
% W_e=[W_e_plus W_e_moins W_e_0];
% Omega_e=inv(W_e);
% 
% Omega_e_plus=Omega_e(1,:);
% Omega_e_moins=Omega_e(2,:);
% 
% 
% Lambda_e_plus= diag(c_e);
% Lambda_e_moins=-Lambda_e_plus;
% 
% B_e=[0 0 1];
% s_e=1;
% 
% M_plus =B_e*M_e*W_e_plus *Lambda_e_plus;
% M_moins=B_e*M_e*W_e_moins*Lambda_e_moins;
% 
% 
% R_e2=-inv(M_moins)*M_plus;
% S_e2= inv(M_moins)*s_e;
% 
% %F_e=A_x_fluid*nx+A_y_fluid*ny;
% F_e2=M_e*(Lambda_e_plus*W_e_plus+Lambda_e_moins*W_e_moins*R_e2)*Omega_e_plus;
% S_e2=M_e*                        Lambda_e_moins*W_e_moins*S_e2;
% 
% 



nx=cos(vec_theta);
ny=sin(vec_theta);
Phi=Phi_fluid_vector(nx,ny,Z_e,Shift_fluid);


II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a,b,[c_edge c_edge]);
MM=kron(II,-F_e);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_edge);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_edge);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;


II=int_edge_1vectorielle(1j*k_e*[nx;ny],a,b,c_edge);
MM=kron(II,S_e);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_edge);  
F(indice_test)=F(indice_test)+Phi.'*MM;

% for i_thetapsi=1:nb_theta
%     theta_psi=vec_theta(i_thetapsi);
%     n_psi=[cos(theta_psi);sin(theta_psi)];
%     Psi_e=conj(Phi_fluid(cos(theta_psi),sin(theta_psi),Z_e));
%     ii=(i_thetapsi-1)+dof_start_element(e_edge)
%     
%     F(ii)=F(ii)+Psi_e'*S_e*...
%         int_edge_1(j*k_e*n_psi,a,b,c_edge);
% end





