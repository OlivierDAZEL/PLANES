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



%! Transmitted fields

for i_T=1:nb_T
    F3=TR6_PW(length_edge,vec_k_x_t(i_T),a);
    
    index_force=dof_A(uy(node));
    index_F_elem=find(index_force);
    index_F_global=index_force(index_F_elem);
    A_2(index_F_global,nb_dof_FEM+nb_R+i_T)=A_2(index_F_global,nb_dof_FEM+nb_R+i_T)-F3(index_F_elem);
    
    A_2(nb_dof_FEM+nb_R+i_T,index_F_global)=A_2(nb_dof_FEM+nb_R+i_T,index_F_global)+F3(index_F_elem)';
end


