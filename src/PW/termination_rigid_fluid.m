% Normal displacement=1
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_2(1))=SV_2(1,1)*exp(-1j*k_z_2*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(2))=SV_2(1,2);