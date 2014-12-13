function Mat_PW=build_global_PW_matrices(k_x,omega,multilayer,termination,nb_layers,nb_amplitudes,n_w,k_air,air)


% Initialization of the matrix
Mat_PW=zeros(nb_amplitudes-1,nb_amplitudes);

% Creation of the equation
number_of_eq=0;
% Space shift for the first layer being x=0
x_interface=-multilayer(1).d;
% Loop on the layers
for i_interface=1:nb_layers-1
    
    %Type of media on both sides and attribution of the dof
    medium_1=multilayer(i_interface).mat;
    dof_medium_1=sum(n_w(1:i_interface-1))+(1:n_w(i_interface));
    medium_2=multilayer(i_interface+1).mat;
    dof_medium_2=sum(n_w(1:i_interface))+(1:n_w(i_interface+1));
    
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
% Last interface
% the last medium is #2 of the end of the loop
if termination==0 % Rigid backing
    switch floor(medium_2/1000)
        case {0 2 3}
            termination_rigid_fluid
        case 1
            termination_rigid_elas
        case 4
            termination_rigid_PEM
    end
else % transmission problem
    % Medium # 1= air last medium is #2
    k_z_1=sqrt(k_air^2-k_x^2);
    SV_1=State_fluid(k_x,k_z_1,air.K);
    switch floor(medium_2/1000)
        case {0 2 3}
            termination_trans_fluid
        case 1
            termination_trans_elas
        case {4 5}
            termination_trans_PEM
    end
    
end

end