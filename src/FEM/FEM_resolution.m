tic
for i_f=1:abs(nb_frequencies)
    FEM_progress=100*i_f/abs(nb_frequencies)
    
    freq=vec_freq(i_f);
    omega=2*pi*freq;
    k_air=omega/air.c;
    create_wave_vectors
    
    % Construction of the linear system
    
    A= sparse(nb_dof_FEM,nb_dof_FEM);
    F=zeros(nb_dof_FEM,1);
    
    
    if (nb_media.acou~=0)
        A=A+H_acou/(air.rho*omega^2)-Q_acou/(air.K);
    end
    
    
    if (nb_media.elas~=0)
        for i_mat=1:nb_media.elas
            eval(['Mat_elas_' num2str(num_media.elas(i_mat))])
            eval(['A=A+(lambda_solide+2*mu_solide)*K0_elas_',num2str(i_mat),'+mu_solide*K1_elas_',num2str(i_mat),'-omega^2*rho_solide*M_elas_',num2str(i_mat),';']);
        end
    end
    
    
    if (nb_media.eqf~=0)
        for i_mat=1:nb_media.eqf
            eval(['Mat_PEM_' num2str(num_media.eqf(i_mat))])
            properties_jca
            eval(['A=A+H_eqf_',num2str(i_mat),'/(rho_eq_til*omega^2)-Q_eqf_',num2str(i_mat),'/(K_eq_til);']);
        end
    end
    
    if (nb_media.limp~=0)
        for i_mat=1:nb_media.limp
            eval(['Mat_PEM_' num2str(num_media.limp(i_mat))])
            properties_limp
            eval(['A=A+H_limp_',num2str(i_mat),'/(rho_limp*omega^2)-Q_limp_',num2str(i_mat),'/(K_eq_til);']);
        end
    end
    if (nb_media.pem98~=0)
        for i_mat=1:nb_media.pem98
            eval(['Mat_PEM_' num2str(num_media.pem98(i_mat))])
            properties_jca
            properties_PEM
            eval(['A=A+P_hat*K0_pem98_',num2str(i_mat),'+N*K1_pem98_',num2str(i_mat),'-omega^2*rho_til*M_pem98_',num2str(i_mat),';']);
            eval(['A=A+H_pem98_',num2str(i_mat),'/(rho_eq_til*omega^2)-Q_pem98_',num2str(i_mat),'/(K_eq_til);']);
            eval(['A=A-gamma_til*(C_pem98_',num2str(i_mat),'+C_pem98_',num2str(i_mat),'.'');']);
        end
    end
    if (nb_media.pem01~=0)
        for i_mat=1:nb_media.pem01
            eval(['Mat_PEM_' num2str(num_media.pem01(i_mat))])
            properties_jca
            properties_PEM
            eval(['A=A+P_hat*K0_pem01_',num2str(i_mat),'+N*K1_pem01_',num2str(i_mat),'-omega^2*rho_til*M_pem01_',num2str(i_mat),';']);
            eval(['A=A+H_pem01_',num2str(i_mat),'/(rho_eq_til*omega^2)-Q_pem01_',num2str(i_mat),'/(K_eq_til);']);
            eval(['A=A-(gamma_til+1)*(C_pem01_',num2str(i_mat),'+C_pem01_',num2str(i_mat),'.'')-(Cp_pem01_',num2str(i_mat),'+Cp_pem01_',num2str(i_mat),'.'');']);
        end
    end
    
    
    if nb_interfaces~=0
        apply_FSI
    end
    
    
    
    if nb_MMT~=0
        
        TT=build_FEM_transfer(k_air*sin(theta_MMT),element_MMT_moins,element_MMT_plus,omega,multilayer_femtmm,k_air,air);
        
        
        switch floor(element_MMT_moins/1000)
            case {0 2 3}
                switch floor(element_MMT_plus/1000)
                    case {0 2 3}
                        link_FEMTMM_fluid_fluid
                end
            case 1
                switch floor(element_MMT_plus/1000)
                    case 1
                        link_FEMTMM_elas_elas
                end
                
                
        end
    end
    
    %     A_JAP=A;
    %     F_JAP=F;
    
    if (nb_R+nb_T)>0
        DtN_application
    end
    
    if length(periodicity)>0
        periodicity_condition_application
    end
    
    
    
    %disp('Resolution of the system')
    X=A\F;
    
    %  X=A_2\F_2;
    
    
    sol(dof_back)=X(1:nb_dof_FEM);
    
    if exist('DtN_plate_R')
        rflx=T_back*  X(nb_dof_FEM+(1:size_info_vector_R*nb_R));
    else
        rflx=X(nb_dof_FEM+(1:size_info_vector_R*nb_R));
    end
    
    if exist('DtN_plate_R')
        trans=T_back_T*X(nb_dof_FEM+size_info_vector_R*nb_R+(1:size_info_vector_T*nb_T));
    else
        trans=X(nb_dof_FEM+size_info_vector_R*nb_R+(1:size_info_vector_T*nb_T));
    end
    
    
    
    if export_profiles==1
        if i_f==1
            sol_export=sol.';
            rflx_export=rflx.';
            trans_export=trans.';
        else
            sol_export=[sol_export sol.'];
            rflx_export=[rflx_export rflx.'];
            trans_export=[trans_export trans.'];
        end
        
        
        
        
    end
    
    
    if export_nrj==1
        Ks(i_f)=(omega^2/4)*rho_1*X(1:nb_dof_FEM)'*M_pem01_1*X(1:nb_dof_FEM);
        Kf(i_f)=(omega^2/4)*(real(rho_f_til)*X(1:nb_dof_FEM)'*M_pem01_1*X(1:nb_dof_FEM)+real(1/conj(rho_eq_til*omega^4))*X(1:nb_dof_FEM)'*H_pem01_1*X(1:nb_dof_FEM)-(2/omega^2)*imag(phi/alpha_til)*imag(X(1:nb_dof_FEM)'*C_pem01_1*X(1:nb_dof_FEM)));
        Wdef(i_f)=(1/4)*(real(P_hat)*X(1:nb_dof_FEM)'*K0_pem01_1*X(1:nb_dof_FEM)+real(N)*X(1:nb_dof_FEM)'*K1_pem01_1*X(1:nb_dof_FEM)+(phi^2*real(R_til)/abs(R_til)^2)*X(1:nb_dof_FEM)'*Q_pem01_1*X(1:nb_dof_FEM));
        
        W_vis(i_f)=(-pi*omega^2)*(imag(rho_til)*X(1:nb_dof_FEM)'*M_pem01_1*X(1:nb_dof_FEM)-imag(1/(rho_eq_til*omega^4))*X(1:nb_dof_FEM)'*H_pem01_1*X(1:nb_dof_FEM)+(2/omega^2)*imag(phi/alpha_til)*real(X(1:nb_dof_FEM)'*C_pem01_1*X(1:nb_dof_FEM)));
        W_struct(i_f)=pi*(imag(P_hat)*X(1:nb_dof_FEM)'*K0_pem01_1*X(1:nb_dof_FEM)+imag(N)*X(1:nb_dof_FEM)'*K1_pem01_1*X(1:nb_dof_FEM));
        W_therm(i_f)=(pi*phi^2*imag(R_til)/abs(R_til)^2)*X(1:nb_dof_FEM)'*Q_pem01_1*X(1:nb_dof_FEM);
        
        W_dis(i_f)=W_vis(i_f)+W_struct(i_f)+W_therm(i_f);
  
        I_inc(i_f)=(period/air.Z)/(2*vec_freq(i_f));
    end
    
    
    R_EF(i_f)=rflx(1);
    abs_EF(i_f)=1-sum(real(vec_k_z).'.*abs(rflx(1:size_info_vector_R:end)).^2)/real(k_z);
     
    %    fprintf(file_abs_id,'%1.15e \t %1.15e \n',freq,abs_EF(i_f));
    
    if nb_T~=0
        TL_EF(i_f)=full(-10*log10(abs(sum(real(vec_k_z_t).'.*abs(trans(1:size_info_vector_T:end)).^2)/real(k_z))));
        fprintf(file_TL_id,'%1.15e \t %1.15e \n',freq,TL_EF(i_f));
    end
    
    
    
    
    
end
time_FEM=toc;

info_FEM


%close (h);


% if export_profiles==1
%     save([name_directory_profiles,'vec_freq'],'vec_freq');
%     save([name_directory_profiles,'elements'],'elements');
%     save([name_directory_profiles,'nodes'],'nodes');
%     save([name_directory_profiles,'element_label'],'element_label');
%     displacements_x=sol_export(1:3:end,:);
%     displacements_y=sol_export(2:3:end,:);
%     pressures=sol_export(3:3:end,:);
%     save([name_directory_profiles,'displacements_x'],'displacements_x');
%     save([name_directory_profiles,'displacements_y'],'displacements_y');
%     save([name_directory_profiles,'pressures'],'pressures');
% end



% x=linspace(0,max(nodes(:,2)),300);
% p_anal=exp(-j*k_air*(x-1-delta_y_MMT))+exp(j*k_air*(x-1-delta_y_MMT));
% p_anal=exp(-j*k_air*x)+exp(j*k_air*x)*exp(-2*j*k_air*(1+delta_y_MMT));
%
% figure
% plot(nodes(:,2),abs(pressures),'.')
% hold on
% plot(x,abs(p_anal),'r')
%
% figure
% plot(nodes(:,2),angle(-pressures),'.')
% hold on
% plot(x,angle(-p_anal),'r')



% figure
% plot(vec_freq,abs_EF,'.')
% hold on
% maine=load('../../Maine/TCLTK/out.dat');
% plot(maine(:,1),maine(:,4))


% figure
% plot(vec_freq,real(R_EF),'.')
% hold on
% plot(vec_freq,imag(R_EF),'r.')
% maine=load('../../Maine/TCLTK/out.dat');
% plot(maine(:,1),maine(:,2))
% plot(maine(:,1),maine(:,3),'r')




% figure
% semilogx(vec_freq,TL,'.')
% hold on
% maine=load('../../Maine/TCLTK/out.dat');
% semilogx(maine(:,1),maine(:,2))

%eval(['rmpath(' list_path ');'])