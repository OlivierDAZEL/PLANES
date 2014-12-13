for ie=1:nb_interfaces
    
    x1=nodes(interfaces(ie,1),1);
    y1=nodes(interfaces(ie,1),2);
    x2=nodes(interfaces(ie,2),1);
    y2=nodes(interfaces(ie,2),2);
    length_edge=sqrt((x2-x1)^2+(y2-y1)^2);
    if (x1<x2)
        a=x1;
        node(1)=interfaces(ie,1);
        node(2)=interfaces(ie,2);
        node(3)=interfaces(ie,6);
    else
        a=x2;
        node(2)=interfaces(ie,1);
        node(1)=interfaces(ie,2);
        node(3)=interfaces(ie,6);
    end
    
    
    a1(1)=nodes(node(1),1);
    a1(2)=nodes(node(1),2);
    a2(1)=nodes(node(2),1);
    a2(2)=nodes(node(2),2);
    
    vec_tangent=a2-a1;
    vec_normal=[vec_tangent(2) vec_tangent(1)];
    vec_normal=vec_normal/norm(vec_normal);
    
    center_e1=mean(nodes(elements(interfaces(ie,3),:),1:2),1);
    center_e2=mean(nodes(elements(interfaces(ie,4),:),1:2),1);
    vec_temp=nodes(interfaces(ie,6),1:2)-center_e1;
    temp=vec_normal*vec_temp';
    if temp<0
        vec_normal=-vec_normal;
    end
    
    if floor(element_label(interfaces(ie,3))/1000)==1
        n_elas=vec_normal;
        n_acou=-vec_normal;
    else
        n_elas=-vec_normal;
        n_acou=vec_normal;
    end
    
    
    FSIe=TR6_FSI(a1,a2);
    
    index_force_p=dof_A(p(node));
    index_F_elem_p=find(index_force_p);
    index_F_global_p=index_force_p(index_F_elem_p);
    
    
    index_force_ux=dof_A(ux(node));
    index_F_elem_ux=find(index_force_ux);
    index_F_global_ux=index_force_ux(index_F_elem_ux);
   
    index_force_uy=dof_A(uy(node));
    index_F_elem_uy=find(index_force_uy);
    index_F_global_uy=index_force_uy(index_F_elem_uy);


    A(index_F_global_p,index_F_global_ux)=A(index_F_global_p,index_F_global_ux)-n_acou(1)*FSIe(index_F_elem_p,index_F_elem_ux);    
    A(index_F_global_p,index_F_global_uy)=A(index_F_global_p,index_F_global_uy)-n_acou(2)*FSIe(index_F_elem_p,index_F_elem_uy);
    A(index_F_global_ux,index_F_global_p)=A(index_F_global_ux,index_F_global_p)-n_elas(1)*FSIe(index_F_elem_ux,index_F_elem_p);    
    A(index_F_global_uy,index_F_global_p)=A(index_F_global_uy,index_F_global_p)-n_elas(2)*FSIe(index_F_elem_uy,index_F_elem_p);    
    
    
end