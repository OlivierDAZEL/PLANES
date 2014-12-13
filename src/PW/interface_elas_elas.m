
% Initialization of the State vectors in media 1 and 2
eval(['Mat_elas_' num2str(medium_1-1000*floor(medium_1/1000))])
delta_P=omega*sqrt(rho_solide/(lambda_solide+2*mu_solide));
delta_s=omega*sqrt(rho_solide/mu_solide);
k_z_1=sqrt([delta_P delta_s].^2-k_x^2);
SV_1=State_elas(k_x,k_z_1,delta_P,delta_s,lambda_solide,mu_solide);

eval(['Mat_elas_' num2str(medium_2-1000*floor(medium_2/1000))])
delta_P=omega*sqrt(rho_solide/(lambda_solide+2*mu_solide));
delta_s=omega*sqrt(rho_solide/mu_solide);
k_z_2=sqrt([delta_P delta_s].^2-k_x^2);
SV_2=State_elas(k_x,k_z_2,delta_P,delta_s,lambda_solide,mu_solide);



% Continuity of all the fields

for i_field=1:4
    number_of_eq=number_of_eq+1;
    Mat_PW(number_of_eq,dof_medium_1(1))= SV_1(i_field,1)*exp(-1j*k_z_1(1)*multilayer(i_interface).d);
    Mat_PW(number_of_eq,dof_medium_1(2))= SV_1(i_field,2)*exp(-1j*k_z_1(2)*multilayer(i_interface).d);
    Mat_PW(number_of_eq,dof_medium_1(3))= SV_1(i_field,3);
    Mat_PW(number_of_eq,dof_medium_1(4))= SV_1(i_field,4);
    
    
    Mat_PW(number_of_eq,dof_medium_2(1))=-SV_2(i_field,1);
    Mat_PW(number_of_eq,dof_medium_2(2))=-SV_2(i_field,2);
    Mat_PW(number_of_eq,dof_medium_2(3))=-SV_2(i_field,3)*exp(-1j*k_z_2(1)*multilayer(i_interface+1).d);
    Mat_PW(number_of_eq,dof_medium_2(4))=-SV_2(i_field,4)*exp(-1j*k_z_2(2)*multilayer(i_interface+1).d);
    
end

