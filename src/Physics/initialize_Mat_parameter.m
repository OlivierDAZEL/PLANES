function Mat_parameter=initialize_Mat_parameter(index_label,index_element,air,omega)

nb_mat=length(index_label);
Mat_parameter=zeros(8,nb_mat);

for ii=1:nb_mat
    typ_elem=floor(index_label(ii)/1000);
    switch typ_elem
        case {0,8}
            Mat_parameter(1,ii)=air.rho;
            Mat_parameter(2,ii)=air.K;
        case {4,5}
            num_mat=index_label(ii)-typ_elem*1000;
            eval(['Mat_porous_' num2str(num_mat)])
            properties_JCA
            properties_PEM
            
            Mat_parameter(1,ii)=rho_eq_til;
            Mat_parameter(2,ii)=K_eq_til;
            Mat_parameter(3,ii)=gamma_til;
            Mat_parameter(4,ii)=A_hat;
            Mat_parameter(5,ii)=P_hat;
            Mat_parameter(6,ii)=N;
            Mat_parameter(7,ii)=rho_s_til;
            Mat_parameter(8,ii)=rho_til;
            Mat_parameter(9,ii)=phi;
        case {2}
            num_mat=index_label(ii)-2000;
            eval(['Mat_porous_' num2str(num_mat)]);
            properties_JCA;
            Mat_parameter(1,ii)=rho_eq_til;
            Mat_parameter(2,ii)=K_eq_til;
        case {1}
            num_mat=index_label(ii)-1000;
            eval(['Mat_elas_' num2str(num_mat)]);
            Mat_parameter(1,ii)=lambda_solide;
            Mat_parameter(2,ii)=mu_solide;
            Mat_parameter(3,ii)=rho_solide;
        case default
            stop666
    end
    
end


