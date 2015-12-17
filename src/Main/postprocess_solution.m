sol(dof_back)=X(1:nb.dof_FEM);

% if ((exist(name.compute_error)==2))
%     eval('eval(name.compute_error)')
% end

if exist('DtN_plate_R')
  %  rflx=T_back*  X(nb.dof_FEM+nb.dof_DGM+(1:size_info_vector_R*nb.R));
    rflx=X(nb.dof_FEM+nb.dof_DGM+(1:size_info_vector_R*nb.R));
else
    rflx=X(nb.dof_FEM+nb.dof_DGM+(1:size_info_vector_R*nb.R));
end

if exist('DtN_plate_R')
   % trans=T_back_T*X(nb.dof_FEM+nb.dof_DGM+size_info_vector_R*nb.R+(1:size_info_vector_T*nb.T));
    trans=X(nb.dof_FEM+nb.dof_DGM+size_info_vector_R*nb.R+(1:size_info_vector_T*nb.T));
else
    trans=X(nb.dof_FEM+nb.dof_DGM+size_info_vector_R*nb.R+(1:size_info_vector_T*nb.T));
end

if export.profiles==1
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


if export.nrj==1
      compute_L2  
%     I_inc(i_f)=(period/air.Z)/(2*frequency.vec(i_f));
%     if num_media.pem01~=0
%         Ks(i_f)=(omega^2/4)*rho_1*X(1:nb.dof_FEM)'*M_pem01_1*X(1:nb.dof_FEM);
%         Kf(i_f)=(omega^2/4)*(real(rho_f_til)*X(1:nb.dof_FEM)'*M_pem01_1*X(1:nb.dof_FEM)+real(1/conj(rho_eq_til*omega^4))*X(1:nb.dof_FEM)'*H_pem01_1*X(1:nb.dof_FEM)-(2/omega^2)*imag(phi/alpha_til)*imag(X(1:nb.dof_FEM)'*C_pem01_1*X(1:nb.dof_FEM)));
%         Wdef(i_f)=(1/4)*(real(P_hat)*X(1:nb.dof_FEM)'*K0_pem01_1*X(1:nb.dof_FEM)+real(N)*X(1:nb.dof_FEM)'*K1_pem01_1*X(1:nb.dof_FEM)+(phi^2*real(R_til)/abs(R_til)^2)*X(1:nb.dof_FEM)'*Q_pem01_1*X(1:nb.dof_FEM));
%         
%         W_vis(i_f)=(-pi*omega^2)*(imag(rho_til)*X(1:nb.dof_FEM)'*M_pem01_1*X(1:nb.dof_FEM)-imag(1/(rho_eq_til*omega^4))*X(1:nb.dof_FEM)'*H_pem01_1*X(1:nb.dof_FEM)+(2/omega^2)*imag(phi/alpha_til)*real(X(1:nb.dof_FEM)'*C_pem01_1*X(1:nb.dof_FEM)));
%         W_struct(i_f)=pi*(imag(P_hat)*X(1:nb.dof_FEM)'*K0_pem01_1*X(1:nb.dof_FEM)+imag(N)*X(1:nb.dof_FEM)'*K1_pem01_1*X(1:nb.dof_FEM));
%         W_therm(i_f)=(pi*phi^2*imag(R_til)/abs(R_til)^2)*X(1:nb.dof_FEM)'*Q_pem01_1*X(1:nb.dof_FEM);
%         W_elas(i_f)=pi*(imag(lambda_solide+2*mu_solide)*X(1:nb.dof_FEM)'*K0_elas_1*X(1:nb.dof_FEM)+imag(mu_solide)*X(1:nb.dof_FEM)'*K1_elas_1*X(1:nb.dof_FEM));
%         
%         abs_vis(i_f)=W_vis(i_f)/I_inc(i_f);
%         abs_struct(i_f)=W_struct(i_f)/I_inc(i_f);
%         abs_therm(i_f)=W_therm(i_f)/I_inc(i_f);
%         abs_elas(i_f)=W_elas(i_f)/I_inc(i_f);
%         
%         abs_dis(i_f)=(abs_vis(i_f)+abs_struct(i_f)+abs_therm(i_f)+abs_elas(i_f));
%         
%     end
end

if nb.R~=0
    reflex_FEM(i_f)=rflx(1);
    abs_EF(i_f)=1-sum(real(vec_k_z).'.*abs(rflx(1:size_info_vector_R:end)).^2)/real(k_z);
end

if nb.T~=0
    TL_EF(i_f)=full(-10*log10(abs(sum(real(vec_k_z_t).'.*abs(trans(1:size_info_vector_T:end)).^2)/real(k_z))));
    fprintf(file_TL_id,'%1.15e \t %1.15e \n',freq,TL_EF(i_f));
end

if profiles.y==1
    plot_sol_PLANES_y
end
if profiles.map==1
        disp('Creation of the 2D plot, this can take time')
        plot_sol_PLANES_map
end
if profiles.custom~=0
    eval(['plot_sol_TR6_custom_' , num2str(profiles.custom)]);
end

% custom plots
if isfield(profiles, 'custom_plots') && length(profiles.custom_plots)~=0
	for i=1:length(profiles.custom_plots)
		eval(profiles.custom_plots{i});
	end
end


clear('sol')
