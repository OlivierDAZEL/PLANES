for ii=1:size(edges_MMT_moins,1)
    
    node_moins=edges_MMT_moins(ii,[1 2 6]);
    node_plus=edges_MMT_plus(ii,[1 2 6]);
    
    
    a1(1)=nodes(node_moins(1),1);
    a1(2)=nodes(node_moins(1),2);
    a2(1)=nodes(node_moins(2),1);
    a2(2)=nodes(node_moins(2),2);
    
    FSIe=TR6_FSI(a1,a2);
    
    
    index_force_p_moins=dof_A(p(node_moins));
    index_F_elem_p_moins=find(index_force_p_moins);
    index_F_global_p_moins=index_force_p_moins(index_F_elem_p_moins);
    
    
    index_force_p_plus=dof_A(p(node_plus));
    index_F_elem_p_plus=find(index_force_p_plus);
    index_F_global_p_plus=index_force_p_plus(index_F_elem_p_plus);
    
    A(index_F_global_p_moins,index_F_global_p_moins)=A(index_F_global_p_moins,index_F_global_p_moins)-TT(1,1)*(FSIe(index_F_elem_p_moins,index_F_elem_p_moins));
    A(index_F_global_p_moins,index_F_global_p_plus) =A(index_F_global_p_moins,index_F_global_p_plus) -TT(1,2)*(FSIe(index_F_elem_p_moins,index_F_elem_p_plus));
    A(index_F_global_p_plus,index_F_global_p_moins) =A(index_F_global_p_plus,index_F_global_p_moins) -TT(2,1)*(FSIe(index_F_elem_p_plus,index_F_elem_p_moins));
    A(index_F_global_p_plus,index_F_global_p_plus)  =A(index_F_global_p_plus,index_F_global_p_plus)  -TT(2,2)*(FSIe(index_F_elem_p_plus,index_F_elem_p_plus));
end