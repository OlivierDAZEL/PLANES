x1=nodes(loads(ie,1),1);
y1=nodes(loads(ie,1),2);
x2=nodes(loads(ie,2),1);
y2=nodes(loads(ie,2),2);
length_edge=sqrt((x2-x1)^2+(y2-y1)^2);
if (x1<x2)
    a=x1;
    node(1)=loads(ie,1);
    node(2)=loads(ie,2);
    node(3)=loads(ie,6);
else
    a=x2;
    node(2)=loads(ie,1);
    node(1)=loads(ie,2);
    node(3)=loads(ie,6);
end
%node
% Incident field
F3=TR6_PW(length_edge,k_x,a);

% Terme p_a delta u_y champs incident
index_force=dof_A(uy(node));
index_F_elem=find(index_force);
index_F_global=index_force(index_F_elem);
F_2(index_F_global)=F_2(index_F_global)+F3(index_F_elem);
%signe ok

% Terme u_a delta p champs incident
index_force=dof_A(p(node));
index_F_elem=find(index_force);
index_F_global=index_force(index_F_elem);
F_2(index_F_global)=F_2(index_F_global)+F3(index_F_elem)*(1i*k_z)/(rho_0*omega^2);
%signe ok

%! Reflected fields

for i_R=1:nb_R
    % Terme p_a delta u_y champs reflechi

    F3=TR6_PW(length_edge,vec_k_x(i_R),a);
    index_force=dof_A(uy(node));
    index_F_elem=find(index_force);
    index_F_global=index_force(index_F_elem);
    A_2(index_F_global,nb_dof_FEM+i_R)=A_2(index_F_global,nb_dof_FEM+i_R)-F3(index_F_elem);
    
    
    % Terme u_a delta p champs reflechi

    index_force=dof_A(p(node));
    index_F_elem=find(index_force);
    index_F_global=index_force(index_F_elem);
    A_2(index_F_global,nb_dof_FEM+i_R)=A_2(index_F_global,nb_dof_FEM+i_R)+F3(index_F_elem)*(1i*vec_k_z(i_R))/(rho_0*omega^2);
    
    %%%%%%% Equation suppl?mentare sur la pression
    A_2(nb_dof_FEM+i_R,index_F_global)=A_2(nb_dof_FEM+i_R,index_F_global)+F3(index_F_elem)';
end



a1(1)=nodes(node(1),1);
a1(2)=nodes(node(1),2);
a2(1)=nodes(node(2),1);
a2(2)=nodes(node(2),2);

FSIe=TR6_FSI(a1,a2);


index_force_p=dof_A(p(node));
index_F_elem_p=find(index_force_p);
index_F_global_p=index_force_p(index_F_elem_p);

index_force_u=dof_A(uy(node));
index_F_elem_u=find(index_force_u);
index_F_global_u=index_force_u(index_F_elem_u);
% Signe ok mercredi
A_2(index_F_global_p,index_F_global_u)=A_2(index_F_global_p,index_F_global_u)-(FSIe(index_F_elem_p,index_F_elem_u));
