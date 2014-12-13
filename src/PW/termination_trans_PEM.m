% Continuity of normal displacement
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_2(1))=-SV_2(3,1)*exp(-1j*k_z_2(1)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(2))=-SV_2(3,2)*exp(-1j*k_z_2(2)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(3))=-SV_2(3,3)*exp(-1j*k_z_2(3)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(4))=-SV_2(3,4);
Mat_PW(number_of_eq,dof_medium_2(5))=-SV_2(3,5);
Mat_PW(number_of_eq,dof_medium_2(6))=-SV_2(3,6);
Mat_PW(number_of_eq,nb_amplitudes)  = SV_1(1,1);
% Continuity of pressure
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_2(1))=-SV_2(5,1)*exp(-1j*k_z_2(1)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(2))=-SV_2(5,2)*exp(-1j*k_z_2(2)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(3))=-SV_2(5,3)*exp(-1j*k_z_2(3)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(4))=-SV_2(5,4);
Mat_PW(number_of_eq,dof_medium_2(5))=-SV_2(5,5);
Mat_PW(number_of_eq,dof_medium_2(6))=-SV_2(5,6);
Mat_PW(number_of_eq,nb_amplitudes)  = SV_1(2,1);

% Nullity of sigma_xx
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_2(1))=SV_2(1,1)*exp(-1j*k_z_2(1)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(2))=SV_2(1,2)*exp(-1j*k_z_2(2)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(3))=SV_2(1,3)*exp(-1j*k_z_2(3)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(4))=SV_2(1,4);
Mat_PW(number_of_eq,dof_medium_2(5))=SV_2(1,5);
Mat_PW(number_of_eq,dof_medium_2(6))=SV_2(1,6);


% Nullity of sigma_zz
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_2(1))=SV_2(4,1)*exp(-1j*k_z_2(1)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(2))=SV_2(4,2)*exp(-1j*k_z_2(2)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(3))=SV_2(4,3)*exp(-1j*k_z_2(3)*multilayer(i_interface+1).d);
Mat_PW(number_of_eq,dof_medium_2(4))=SV_2(4,4);
Mat_PW(number_of_eq,dof_medium_2(5))=SV_2(4,5);
Mat_PW(number_of_eq,dof_medium_2(6))=SV_2(4,6);