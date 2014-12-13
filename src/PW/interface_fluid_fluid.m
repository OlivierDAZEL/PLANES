% Initialization of the State vectors in media 1 and 2
if medium_1==0
    k_z_1=sqrt(k_air^2-k_x^2);
    SV_1=State_fluid(k_x,k_z_1,air.K);
else
    eval(['Mat_fluid_' num2str(medium_1-1000*floor(medium_1/1000))]);
    properties_jca
    k_z_1=sqrt((omega*sqrt(rho_eq_til/K_eq_til))^2-k_x^2);
    SV_1=State_fluid(k_x,k_z_1,K_eq_til);
end
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
Mat_PW(number_of_eq,dof_medium_1(1))= SV_1(1,1)*exp(-1j*k_z_1*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(2))= SV_1(1,2);
Mat_PW(number_of_eq,dof_medium_2(1))=-SV_2(1,1);
Mat_PW(number_of_eq,dof_medium_2(2))=-SV_2(1,2)*exp(-1j*k_z_2*multilayer(i_interface+1).d);
% Continuity of the pressure
number_of_eq=number_of_eq+1;
Mat_PW(number_of_eq,dof_medium_1(1))= SV_1(2,1)*exp(-1j*k_z_1*multilayer(i_interface).d);
Mat_PW(number_of_eq,dof_medium_1(2))= SV_1(2,2);
Mat_PW(number_of_eq,dof_medium_2(1))=-SV_2(2,1);
Mat_PW(number_of_eq,dof_medium_2(2))=-SV_2(2,2)*exp(-1j*k_z_2*multilayer(i_interface+1).d);