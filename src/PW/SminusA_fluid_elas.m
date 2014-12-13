% Initialization of the State vectors in media 1 and 2

eval(['Mat_elas_' num2str(multilayer(1).mat-1000*floor(multilayer(1).mat/1000))])
delta_P=omega*sqrt(rho_solide/(lambda_solide+2*mu_solide));
delta_s=omega*sqrt(rho_solide/mu_solide);
k_z_2=sqrt([delta_P delta_s].^2-k_x^2);
SV_2=State_elas(k_x,k_z_2,delta_P,delta_s,lambda_solide,mu_solide);
dof_medium_2=nS_minus+[1:4];

% Continuity of normal displacement
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,1)= -1;
Mat_PW(number_of_eq,dof_medium_2(1))=SV_2(2,1);
Mat_PW(number_of_eq,dof_medium_2(2))=SV_2(2,2);
Mat_PW(number_of_eq,dof_medium_2(3))=SV_2(2,3)*exp(-1j*k_z_2(1)*multilayer(1).d);
Mat_PW(number_of_eq,dof_medium_2(4))=SV_2(2,4)*exp(-1j*k_z_2(2)*multilayer(1).d);
% Continuity of the pressure
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,2)= -1;
Mat_PW(number_of_eq,dof_medium_2(1))=-SV_2(3,1);
Mat_PW(number_of_eq,dof_medium_2(2))=-SV_2(3,2);
Mat_PW(number_of_eq,dof_medium_2(3))=-SV_2(3,3)*exp(-1j*k_z_2(1)*multilayer(1).d);
Mat_PW(number_of_eq,dof_medium_2(4))=-SV_2(3,4)*exp(-1j*k_z_2(2)*multilayer(1).d);

% Nullity of sigma_xz
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_2(1))=SV_2(1,1);
Mat_PW(number_of_eq,dof_medium_2(2))=SV_2(1,2);
Mat_PW(number_of_eq,dof_medium_2(3))=SV_2(1,3)*exp(-1j*k_z_2(1)*multilayer(1).d);
Mat_PW(number_of_eq,dof_medium_2(4))=SV_2(1,4)*exp(-1j*k_z_2(2)*multilayer(1).d);



