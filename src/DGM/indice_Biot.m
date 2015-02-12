function ind=indice_Biot_2D(i_elem,i_theta,i_onde,dof_start_element)




ind=i_onde+(i_theta-1)*3+dof_start_element(i_elem)-1;


% ondes dans angles dans elements