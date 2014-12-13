function T=build_FEM_transfer(k_x,element_MMT_minus,element_MMT_plus,omega,multilayer,k_air,air)


nb_layers=length(multilayer);
termination=0;
compute_number_PW_TMM



switch floor(element_MMT_minus/1000)
    
    case 0
        k_z_minus=sqrt(k_air^2-k_x^2);
        SV_minus=State_fluid(k_x,k_z_minus,air.K);
        nS_minus=2;
        dof_FEM=[2];
        boundary_FEM=[1];
    case 1
        eval(['Mat_elas_' num2str(element_MMT_minus-1000*floor(element_MMT_minus/1000))])
        delta_P=omega*sqrt(rho_solide/(lambda_solide+2*mu_solide));
        delta_s=omega*sqrt(rho_solide/mu_solide);
        k_z_minus=sqrt([delta_P delta_s].^2-k_x^2);
        SV_minus=State_elas(k_x,k_z_minus,delta_P,delta_s,lambda_solide,mu_solide);    
        nS_minus=4;
        dof_FEM=[4 2];
        boundary_FEM=[1 3];
    case 2
        eval(['Mat_fluid_' num2str(element_MMT_minus-1000*floor(element_MMT_minus/1000))]);
        properties_jca
        k_z_minus=sqrt((omega*sqrt(rho_eq_til/K_eq_til))^2-k_x^2);
        SV_minus=State_fluid(k_x,k_z_minus,K_eq_til);
        nS_minus=2;
        dof_FEM=[2];
        boundary_FEM=[1];
    case {4 5}
        eval(['Mat_PEM_' num2str(element_MMT_minus-1000*floor(element_MMT_minus/1000))]);
        properties_jca
        properties_PEM
        compute_Biot_waves
        SV_minus=State_PEM(k_x,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til);
        nS_minus=6;
end

switch floor(element_MMT_plus/1000)
    
    case 0
        k_z_plus=sqrt(k_air^2-k_x^2);
        SV_plus=State_fluid(k_x,k_z_plus,air.K);
        nS_plus=2;
        dof_FEM=[dof_FEM nS_minus+[2]];
        boundary_FEM=[boundary_FEM nS_minus+[1]];
        normale_plus=diag([-1]);
    case 1
         eval(['Mat_elas_' num2str(element_MMT_plus-1000*floor(element_MMT_plus/1000))])
        delta_P=omega*sqrt(rho_solide/(lambda_solide+2*mu_solide));
        delta_s=omega*sqrt(rho_solide/mu_solide);
        k_z_plus=sqrt([delta_P delta_s].^2-k_x^2);
        SV_plus=State_elas(k_x,k_z_plus,delta_P,delta_s,lambda_solide,mu_solide);    
        nS_plus=4;
        dof_FEM=[dof_FEM nS_minus+[4 2]];
        boundary_FEM=[boundary_FEM nS_minus+[1 3]];
        normale_plus=diag([-1 -1]);
    case 2
        eval(['Mat_fluid_' num2str(element_MMT_plus-1000*floor(element_MMT_plus/1000))]);
        properties_jca
        k_z_plus=sqrt((omega*sqrt(rho_eq_til/K_eq_til))^2-k_x^2);
        SV_plus=State_fluid(k_x,k_z_plus,K_eq_til);
        nS_plus=2;
        dof_FEM=[dof_FEM nS_minus+[2]];
        boundary_FEM=[boundary_FEM nS_minus+[1]];
        normale_plus=diag([-1]);
    case {4 5}
        eval(['Mat_PEM_' num2str(element_MMT_plus-1000*floor(element_MMT_plus/1000))]);
        properties_jca
        properties_PEM
        compute_Biot_waves
        SV_plus=State_PEM(k_x,k_z_plus,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til);
        nS_plus=6;
end


% Initialization of the matrix
Mat_PW=zeros(nS_minus/2+nb_amplitudes+nS_plus/2,nS_minus+nb_amplitudes+nS_plus);
% Creation of the equation

number_of_eq=0;

if ismember(floor(element_MMT_minus/1000),[0 2 3])
    switch floor(multilayer(1).mat/1000)
        case {0 2 3}
            SminusA_fluid_fluid
        case 1
            SminusA_fluid_elas
        case {4 5}
            SminusA_fluid_PEM
    end
elseif floor(element_MMT_minus/1000)==1
    switch floor(multilayer(1).mat/1000)
        case {0 2 3}
            fhgfghhghgfgfh
        case 1
            fhgfghhghgfgfh
        case {4 5}
            SminusA_PEM_elas
    end
elseif ismember(floor(element_MMT_minus/1000),[4 5])
    switch floor(multilayer(1).mat/1000)
        case {0 2 3}
            fhgfghhghgfgfh
        case 1
            fhgfghhghgfgfh
        case {4 5}
            fhgfghhghgfgfh
    end
end


if ismember(floor(element_MMT_plus/1000),[0 2 3])
    switch floor(multilayer(end).mat/1000)
        case {0 2 3}
            ASplus_fluid_fluid
        case 1
            ASplus_elas_fluid
        case {4 5}
            ASplus_PEM_fluid
    end
elseif floor(element_MMT_plus/1000)==1
    switch floor(multilayer(end).mat/1000)
        case {0 2 3}
            fhgfghhghgfgfh
        case 1
            fhgfghhghgfgfh
        case {4 5}
            ASplus_PEM_elas
    end
elseif ismember(floor(element_MMT_plus/1000),[4 5])
    switch floor(multilayer(end).mat/1000)
        case {0 2 3}
            fhgfghhghgfgfh
        case 1
            fhgfghhghgfgfh
        case {4 5}
            fhgfghhghgfgfh
    end
end



for i_interface=1:nb_layers-1
    
    %Type of media on both sides and attribution of the dof
    medium_1=multilayer(i_interface).mat;
    dof_medium_1=nS_minus+sum(n_w(1:i_interface-1))+(1:n_w(i_interface));
    medium_2=multilayer(i_interface+1).mat;
    dof_medium_2=nS_minus+sum(n_w(1:i_interface))+(1:n_w(i_interface+1));
    
    if ismember(floor(medium_1/1000),[0 2 3])
        switch floor(medium_2/1000)
            case {0 2 3}
                interface_fluid_fluid
            case 1
                interface_fluid_elas
            case {4 5}
                interface_fluid_PEM
        end
    elseif floor(medium_1/1000)==1
        switch floor(medium_2/1000)
            case {0 2 3}
                interface_elas_fluid
            case 1
                interface_elas_elas
            case {4 5}
                interface_elas_PEM
        end
    elseif ismember(floor(medium_1/1000),[4 5])
        switch floor(medium_2/1000)
            case {0 2 3}
                interface_PEM_fluid
            case 1
                interface_PEM_elas
            case {4 5}
                interface_PEM_PEM
        end
    end
end



S_moins=1:nS_minus;
dof_amplitudes=nS_minus+[1:nb_amplitudes];
S_plus=nS_minus+nb_amplitudes+[1:nS_plus];


M11=Mat_PW([1:nS_minus/2],S_moins);
M12=Mat_PW([1:nS_minus/2],dof_amplitudes);
M13=Mat_PW([1:nS_minus/2],S_plus);
M21=Mat_PW(nS_minus/2+[1:nb_amplitudes],S_moins);
M22=Mat_PW(nS_minus/2+[1:nb_amplitudes],dof_amplitudes);
M23=Mat_PW(nS_minus/2+[1:nb_amplitudes],S_plus);
M31=Mat_PW(nS_minus/2+nb_amplitudes+[1:nS_plus/2],S_moins);
M32=Mat_PW(nS_minus/2+nb_amplitudes+[1:nS_plus/2],dof_amplitudes);
M33=Mat_PW(nS_minus/2+nb_amplitudes+[1:nS_plus/2],S_plus);

Gamma=-inv(M22)*[M21 M23];

GGamma=[[M11 M13]+M12*Gamma;[M31 M33]+M32*Gamma];

M_b=GGamma(:,boundary_FEM);
M_d=GGamma(:,dof_FEM);

T=-inv(M_b)*M_d;

T(nS_minus/2+[1:nS_plus/2],:)=normale_plus*T(nS_minus/2+[1:nS_plus/2],:);


