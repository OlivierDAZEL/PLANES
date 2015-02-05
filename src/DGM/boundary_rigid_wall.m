%%%%% Coordinates of the edge's vertices
coord_edge(1:2,1)=nodes(dirichlets(ie,1),1:2)';
coord_edge(1:2,2)=nodes(dirichlets(ie,2),1:2)';

a=coord_edge(:,1);
b=coord_edge(:,2);

h=norm(b-a);
n_ab=(b-a)/h;

%%%%% Element linked to the edge

e_edge=dirichlets(ie,3);
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


W_e_plus=Phi_fluid(nx,ny,Z_e);
W_e_moins=Phi_fluid(-nx,-ny,Z_e);
W_e_0=Phi_fluid_0(nx,ny);


W_e=[W_e_plus W_e_moins W_e_0];
Omega_e=inv(W_e);

Omega_e_plus=Omega_e(1,:);
Omega_e_moins=Omega_e(2,:);


Lambda_e_plus=diag(c_e);
Lambda_e_moins=-Lambda_e_plus;

B_e=[0 0 1];

M_plus =B_e*M_e*W_e_plus *Lambda_e_plus;
M_moins=B_e*M_e*W_e_moins*Lambda_e_moins;

R_e=-inv(M_moins)*M_plus;

%F_e=A_x_fluid*nx+A_y_fluid*ny;
F_e=M_e*(Lambda_e_plus*W_e_plus+Lambda_e_moins*W_e_moins*R_e)*Omega_e_plus;

delta_test=k_e;
delta_champs=k_e;

nx=cos(vec_theta);
ny=sin(vec_theta);
Phi=Phi_fluid_vector(nx,ny,air.Z,Shift_fluid);

II=int_edge_2vectorielle(1j*delta_test*[nx;ny],-1j*delta_champs*[nx;ny],a,b,[c_edge c_edge]);
MM=kron(II,F_e);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_edge);  
indice_champs=((1:nb_theta)-1)+dof_start_element(e_edge);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;

