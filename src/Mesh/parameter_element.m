
if element_label(e_edge)==0
    c_e=c_0;
    k_e=omega/c_e;
    Z_e=Z_0;
    tau_x=1;
    tau_y=1;
    M_e=diag([rho_0,rho_0,1/K_0]);
    
elseif floor(element_label(e_edge)/1000)==2
    
    temp=load('../inputs/materials/porous/1.biot');
    
    
    phi=temp(1);
	sig=temp(2);
	alpha=temp(3);
	LCV=temp(4);
	LCT=temp(5);
	young=temp(6);
	poisson=temp(7);
    eta=temp(8);
	rho_1=temp(9);
    
    tau_x=1;
    tau_y=1;
    calcul_propriete
    calcul_theorique
    
    k_e=delta_eq;
    c_e=omega/delta_eq;
    Z_e=rho_eq_til*c_e;
    M_e=diag([rho_eq_til,rho_eq_til,1/K_eq_til]); 
    

elseif floor(element_label(e_edge)/1000)==8    
    c_e=c_0;
    
    tau_x=1;
    tau_y=1;
    temp=(element_label(e_edge)-8000);
    if (floor(temp/100)==1)
        tau_x=exp(j*pi/4);
    end
    temp=temp-100*floor(temp/100);
    if (floor(temp/10)==1)
        tau_y=exp(j*pi/4);
    end
    

    k_e=omega/c_e;
    Z_e=Z_0;
    
    M_e=diag([rho_0,rho_0,1/K_0]);    
    
elseif floor(element_label(e_edge)/1000)==4
    
    temp=load('../inputs/materials/porous/1.biot');
    tau_x=1;
    tau_y=1;
    
    phi=temp(1);
	sig=temp(2);
	alpha=temp(3);
	LCV=temp(4);
	LCT=temp(5);
	young=temp(6);
	poisson=temp(7);
    eta=temp(8);
	rho_1=temp(9);
    calcul_propriete
    calcul_theorique
    
    M_e=zeros(8,8);
    M_e=[rho_s_til 0 gamma_til*rho_eq_til 0; 0 rho_s_til 0 gamma_til*rho_eq_til;gamma_til*rho_eq_til 0 rho_eq_til 0 ;0 gamma_til*rho_eq_til 0 rho_eq_til];
    M_e(5,5)=1/(A_hat+N);
    M_e(6,6)=1/N;
    M_e(7,7)=1/N;
    M_e(8,8)=1/K_eq_til;
   
else
    disp('Subroutine parameter_element')
    disp('Unknwon fluid type of element')
    stop
end