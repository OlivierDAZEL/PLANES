
% Initialization of the State vectors in media 1 and 2

eval(['Mat_PEM_' num2str(medium_1-1000*floor(medium_1/1000))])
properties_jca
properties_PEM
compute_Biot_waves
k_z_1=sqrt([delta_1 delta_2 delta_3].^2-k_x^2);
SV_1=State_PEM(k_x,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til);


medium_2

if medium_2==0
    k_z_2=sqrt(k_air^2-k_x^2);
    SV_2=State_fluid(k_x,k_z_2,air.K);
else
    eval(['Mat_fluid_' num2str(medium_2-1000*floor(medium_2/1000))])
    properties_jca
    k_z_2=sqrt((omega*sqrt(rho_eq_til/K_eq_til))^2-k_x^2);
    SV_2=State_fluid(k_x,k_z_2,K_eq_til);
end




% Continuity of normal displacement

number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_1(1))=-SV_1(3,1)*exp(-1j*k_z_1(1)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(2))=-SV_1(3,2)*exp(-1j*k_z_1(2)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(3))=-SV_1(3,3)*exp(-1j*k_z_1(3)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(4))=-SV_1(3,4);
Mat_PW(number_of_eq,dof_medium_1(5))=-SV_1(3,5);
Mat_PW(number_of_eq,dof_medium_1(6))=-SV_1(3,6);
Mat_PW(number_of_eq,dof_medium_2(1))= SV_2(1,1);
Mat_PW(number_of_eq,dof_medium_2(2))= SV_2(1,2)*exp(-1j*k_z_2*multilayer(i_interface+1).d);



% Continuity of pressure
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_1(1))=-SV_1(5,1)*exp(-1j*k_z_1(1)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(2))=-SV_1(5,2)*exp(-1j*k_z_1(2)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(3))=-SV_1(5,3)*exp(-1j*k_z_1(3)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(4))=-SV_1(5,4);
Mat_PW(number_of_eq,dof_medium_1(5))=-SV_1(5,5);
Mat_PW(number_of_eq,dof_medium_1(6))=-SV_1(5,6);
Mat_PW(number_of_eq,dof_medium_2(1))= SV_2(2,1);
Mat_PW(number_of_eq,dof_medium_2(2))= SV_2(2,2)*exp(-1j*k_z_2*multilayer(i_interface+1).d);


% Nullity of sigma_xx
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_1(1))=SV_1(1,1)*exp(-1j*k_z_1(1)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(2))=SV_1(1,2)*exp(-1j*k_z_1(2)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(3))=SV_1(1,3)*exp(-1j*k_z_1(3)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(4))=SV_1(1,4);
Mat_PW(number_of_eq,dof_medium_1(5))=SV_1(1,5);
Mat_PW(number_of_eq,dof_medium_1(6))=SV_1(1,6);


% Nullity of sigma_zz
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_1(1))=SV_1(4,1)*exp(-1j*k_z_1(1)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(2))=SV_1(4,2)*exp(-1j*k_z_1(2)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(3))=SV_1(4,3)*exp(-1j*k_z_1(3)*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(4))=SV_1(4,4);
Mat_PW(number_of_eq,dof_medium_1(5))=SV_1(4,5);
Mat_PW(number_of_eq,dof_medium_1(6))=SV_1(4,6);


