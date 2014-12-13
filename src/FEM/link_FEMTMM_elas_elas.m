for ii=1:size(edges_MMT_moins,1)
    
 
    
    
    node_moins=edges_MMT_moins(ii,[1 2 6]);
    node_plus=edges_MMT_plus(ii,[1 2 6]);
        
    a1(1)=nodes(node_moins(1),1);
    a1(2)=nodes(node_moins(1),2);
    a2(1)=nodes(node_moins(2),1);
    a2(2)=nodes(node_moins(2),2);
    
    
    vec_tangent=a2-a1;
    vec_normal=[vec_tangent(2) vec_tangent(1)];
    vec_normal=vec_normal/norm(vec_normal);
    
    center_element=mean(nodes(elements(edges_MMT_moins(ii,3),:),1:2),1);
    vec_temp=nodes(node_moins(3),1:2)-center_element;
    temp=vec_normal*vec_temp';
    if temp<0
        vec_normal=-vec_normal;
    end
    theta_rot=angle(vec_normal(1)+1i*vec_normal(2))-pi/2  ;
    Mat_rot=[cos(theta_rot) sin(theta_rot);-sin(theta_rot) cos(theta_rot)];
    Mat_rot=[Mat_rot 0*Mat_rot;0*Mat_rot Mat_rot];
    Mat_rotm1=[cos(theta_rot) -sin(theta_rot);sin(theta_rot) cos(theta_rot)];
    Mat_rotm1=[Mat_rotm1 0*Mat_rotm1;0*Mat_rotm1 Mat_rotm1];
    TTrot=Mat_rotm1*TT*Mat_rot;
    
    
    FSIe=TR6_FSI(a1,a2);
    
    index_force_ux_moins=dof_A(ux(node_moins));
    index_F_elem_ux_moins=find(index_force_ux_moins);
    index_F_global_ux_moins=index_force_ux_moins(index_F_elem_ux_moins);   
    
    index_force_uy_moins=dof_A(uy(node_moins));
    index_F_elem_uy_moins=find(index_force_uy_moins);
    index_F_global_uy_moins=index_force_uy_moins(index_F_elem_uy_moins);
    
    index_force_ux_plus=dof_A(ux(node_plus));
    index_F_elem_ux_plus=find(index_force_ux_plus);
    index_F_global_ux_plus=index_force_ux_plus(index_F_elem_ux_plus);
    
    
    index_force_uy_plus=dof_A(uy(node_plus));
    index_F_elem_uy_plus=find(index_force_uy_plus);
    index_F_global_uy_plus=index_force_uy_plus(index_F_elem_uy_plus);
    
    
    
    A(index_F_global_ux_moins,index_F_global_ux_moins)=A(index_F_global_ux_moins,index_F_global_ux_moins)-TTrot(1,1)*(FSIe(index_F_elem_ux_moins,index_F_elem_ux_moins));
    A(index_F_global_ux_moins,index_F_global_uy_moins)=A(index_F_global_ux_moins,index_F_global_uy_moins)-TTrot(1,2)*(FSIe(index_F_elem_ux_moins,index_F_elem_uy_moins));
    A(index_F_global_ux_moins,index_F_global_ux_plus )=A(index_F_global_ux_moins,index_F_global_ux_plus )-TTrot(1,3)*(FSIe(index_F_elem_ux_moins,index_F_elem_ux_plus ));
    A(index_F_global_ux_moins,index_F_global_uy_plus )=A(index_F_global_ux_moins,index_F_global_uy_plus )-TTrot(1,4)*(FSIe(index_F_elem_ux_moins,index_F_elem_uy_plus ));
    
    A(index_F_global_uy_moins,index_F_global_ux_moins)=A(index_F_global_uy_moins,index_F_global_ux_moins)-TTrot(2,1)*(FSIe(index_F_elem_uy_moins,index_F_elem_ux_moins));
    A(index_F_global_uy_moins,index_F_global_uy_moins)=A(index_F_global_uy_moins,index_F_global_uy_moins)-TTrot(2,2)*(FSIe(index_F_elem_uy_moins,index_F_elem_uy_moins));
    A(index_F_global_uy_moins,index_F_global_ux_plus )=A(index_F_global_uy_moins,index_F_global_ux_plus )-TTrot(2,3)*(FSIe(index_F_elem_uy_moins,index_F_elem_ux_plus ));
    A(index_F_global_uy_moins,index_F_global_uy_plus )=A(index_F_global_uy_moins,index_F_global_uy_plus )-TTrot(2,4)*(FSIe(index_F_elem_uy_moins,index_F_elem_uy_plus ));

    A(index_F_global_ux_plus ,index_F_global_ux_moins)=A(index_F_global_ux_plus ,index_F_global_ux_moins)-TTrot(3,1)*(FSIe(index_F_elem_ux_plus ,index_F_elem_ux_moins));
    A(index_F_global_ux_plus ,index_F_global_uy_moins)=A(index_F_global_ux_plus ,index_F_global_uy_moins)-TTrot(3,2)*(FSIe(index_F_elem_ux_plus ,index_F_elem_uy_moins));
    A(index_F_global_ux_plus ,index_F_global_ux_plus )=A(index_F_global_ux_plus ,index_F_global_ux_plus )-TTrot(3,3)*(FSIe(index_F_elem_ux_plus ,index_F_elem_ux_plus ));
    A(index_F_global_ux_plus ,index_F_global_uy_plus )=A(index_F_global_ux_plus ,index_F_global_uy_plus )-TTrot(3,4)*(FSIe(index_F_elem_ux_plus ,index_F_elem_uy_plus ));
    
    A(index_F_global_uy_plus ,index_F_global_ux_moins)=A(index_F_global_uy_plus ,index_F_global_ux_moins)-TTrot(4,1)*(FSIe(index_F_elem_uy_plus ,index_F_elem_ux_moins));
    A(index_F_global_uy_plus ,index_F_global_uy_moins)=A(index_F_global_uy_plus ,index_F_global_uy_moins)-TTrot(4,2)*(FSIe(index_F_elem_uy_plus ,index_F_elem_uy_moins));
    A(index_F_global_uy_plus ,index_F_global_ux_plus )=A(index_F_global_uy_plus ,index_F_global_ux_plus )-TTrot(4,3)*(FSIe(index_F_elem_uy_plus ,index_F_elem_ux_plus ));
    A(index_F_global_uy_plus ,index_F_global_uy_plus )=A(index_F_global_uy_plus ,index_F_global_uy_plus )-TTrot(4,4)*(FSIe(index_F_elem_uy_plus ,index_F_elem_uy_plus ));
   

    
    end