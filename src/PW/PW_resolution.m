for i_f=1:abs(nb_frequencies)
    omega=2*pi*vec_freq(i_f);
    k_air=omega/air.c;
    k_x=k_air*sin(theta);
    
    Mat_PW=build_global_PW_matrices(k_x,omega,multilayer,termination,nb_layers,nb_amplitudes,n_w,k_air,air);
    
    F_PW=-Mat_PW(:,1);
    Mat_PW(:,1)=[];
    X_PW=Mat_PW\F_PW;
    
    abs_PW(i_f)=1-abs(X_PW(1))^2;
    if termination~=0
        TL_PW(i_f)=-20*log10(abs(X_PW(end)));
    end
    
end