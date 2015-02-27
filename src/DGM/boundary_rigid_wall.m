%%%%% Coordinates of the edge's vertices
coord_edge(1:2,1)=nodes(dirichlets(ie,1),1:2)';
coord_edge(1:2,2)=nodes(dirichlets(ie,2),1:2)';

a=coord_edge(:,1);
b=coord_edge(:,2);

%%%%% Element linked to the edge

e_edge=dirichlets(ie,3);
c_edge=mean(nodes(elements(e_edge,:),1:2))';

parameter_element

%%%%% vector normal pointing outwards \Omega_e

centre_edge=(a+b)/2;
n_centre=c_edge-centre_edge;
ne=normal_edge(coord_edge);
if (n_centre'*ne>0)
    ne=-ne;
end
nx=ne(1);
ny=ne(2);


F_e=zeros(3,3);
F_e(1,1)=Z_e*(nx^2);
F_e(1,2)=Z_e*(nx*ny);
F_e(1,3)=-nx;
F_e(2,1)=Z_e*(nx*ny);
F_e(2,2)=Z_e*(ny^2);
F_e(2,3)=-ny;

nx=cos(vec_theta);
ny=sin(vec_theta);
Phi=Phi_fluid_vector(nx,ny,Z_e,Shift_fluid);

II=int_edge_2vectorielle(1j*k_e*[nx;ny],-1j*k_e*[nx;ny],a,b,[c_edge c_edge]);
MM=kron(II,F_e);
indice_test  =((1:nb_theta)-1)+dof_start_element(e_edge);
indice_champs=((1:nb_theta)-1)+dof_start_element(e_edge);
A(indice_test,indice_champs)=A(indice_test,indice_champs)+Phi.'*MM*Phi;
