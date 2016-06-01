% build_FEM_transfer_2D.m
%
% Copyright (C) 2014 < Olivier DAZEL <olivier.dazel@univ-lemans.fr> >
%
% This file is part of PLANES.
%
% PLANES (Porous LAum NumErical Simulator) is a software to compute the
% vibroacoustic response of sound packages containing coupled
% acoustic/elastic/porous substructures. It is mainly based on the
% Finite-Element Method and some numerical methods developped at
% LAUM (http://laum.univ-lemans.fr).
%
% You can download the latest version of PLANES at
% https://github.com/OlivierDAZEL/PLANES
% or find more details on Olivier's webpage
% http://perso.univ-lemans.fr/~odazel/
%
% For any questions or if you want to
% contribute to this project, contact Olivier.
%
% PLANES is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%%


function [ZOD_mat,nb_dof_minus,nb_dof_plus]=build_FEM_transfer_2D(k_x,element_MMT_minus,element_MMT_plus,omega,multilayer,k_air,air)

nb_layers=length(multilayer);
multilayer(1,1).nb=1;
nb_multilayers=1;
multilayer(1,1).termination=0;
compute_number_PW_2D


switch floor(element_MMT_minus/1000)
    
    case 0
        k_z_minus=sqrt(k_air^2-k_x^2);
        SV_minus=State_fluid(k_x,k_z_minus,air.K);
        nS_minus=2;
        dof_FEM=[2];
        boundary_FEM=[1];
        index_ZOD=[3];
    case 1
        eval(['Mat_elas_' num2str(element_MMT_minus-1000*floor(element_MMT_minus/1000))])
        delta_P=omega*sqrt(rho_solide/(lambda_solide+2*mu_solide));
        delta_s=omega*sqrt(rho_solide/mu_solide);
        
        k_z_minus=sqrt([delta_P delta_s].^2-k_x^2);
        SV_minus=State_elas_2D(k_x,delta_P,delta_s,lambda_solide,mu_solide);
        nS_minus=4;
        dof_FEM=[4 2];
        boundary_FEM=[1 3];
        index_ZOD=[1;2];
    case {2 3}
        eval(['Mat_porous_' num2str(element_MMT_minus-1000*floor(element_MMT_minus/1000))]);
        properties_eqf
        k_z_minus=sqrt((omega*sqrt(rho_eq_til/K_eq_til))^2-k_x^2);
        nS_minus=2;
        dof_FEM=[2];
        boundary_FEM=[1];
        index_ZOD=[3];
    case {4 5}
        eval(['Mat_porous_' num2str(element_MMT_minus-1000*floor(element_MMT_minus/1000))]);
        properties_eqf
        properties_PEM
        compute_Biot_waves
        SV_minus=State_PEM_2D(k_x,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til);
        nS_minus=6;
        dof_FEM=[6 2 5];
        boundary_FEM=[1 4 3];
        index_ZOD=[1;2;3];
end

switch floor(element_MMT_plus/1000)
    
    case 0
        k_z_plus=sqrt(k_air^2-k_x^2);
        SV_plus=State_fluid_2D(k_x,k_z_plus,air.K);
        nS_plus=2;
        dof_FEM=[dof_FEM nS_minus+[2]];
        boundary_FEM=[boundary_FEM nS_minus+[1]];
        normale_plus=diag([-1]);
        index_ZOD=[index_ZOD; 3+[3]];
    case 1
        eval(['Mat_elas_' num2str(element_MMT_plus-1000*floor(element_MMT_plus/1000))])
        delta_P=omega*sqrt(rho_solide/(lambda_solide+2*mu_solide));
        delta_s=omega*sqrt(rho_solide/mu_solide);
        SV_plus=State_elas_2D(k_x,delta_P,delta_s,lambda_solide,mu_solide);
        nS_plus=4;
        dof_FEM=[dof_FEM nS_minus+[4 2]];
        boundary_FEM=[boundary_FEM nS_minus+[1 3]];
        normale_plus=diag([-1 -1]);
        index_ZOD=[index_ZOD; 3+[1;2]];
    case {2 3}
        eval(['Mat_porous_' num2str(element_MMT_plus-1000*floor(element_MMT_plus/1000))]);
        properties_eqf
        k_z_plus=sqrt((omega*sqrt(rho_eq_til/K_eq_til))^2-k_x^2);
        nS_plus=2;
        dof_FEM=[dof_FEM nS_minus+[2]];
        boundary_FEM=[boundary_FEM nS_minus+[1]];
        normale_plus=diag([-1]);
        index_ZOD=[index_ZOD; 3+[3]];
    case {4 5}
        eval(['Mat_porous_' num2str(element_MMT_plus-1000*floor(element_MMT_plus/1000))]);
        properties_eqf
        properties_PEM
        compute_Biot_waves
        SV_plus=State_PEM_2D(k_x,k_z_plus,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til);
        nS_plus=6;
        dof_FEM=[dof_FEM nS_minus+[6 2 5]];
        boundary_FEM=[boundary_FEM nS_minus+[1 4 3]];
        normale_plus=diag([-1 -1 -1]);
        index_ZOD=[index_ZOD; 3+[1;2;3]];
end


% Initialization of the matrix
Mat_PW=zeros(nS_minus/2+nb_amplitudes+nS_plus/2,nS_minus+nb_amplitudes+nS_plus);
% Creation of the equation

number_of_eq=0;

if ismember(floor(element_MMT_minus/1000),[0 2 3])
    switch floor(multilayer(1).mat/1000)
        case {0 2 3}
            SminusA_fluid_fluid_2D
        case 1
            SminusA_fluid_elas_2D
        case {4 5}
            SminusA_fluid_PEM_2D
    end
elseif floor(element_MMT_minus/1000)==1
    switch floor(multilayer(1).mat/1000)
        case {0 2 3}
            fhgfghhghgfgfh
        case 1
            SminusA_elas_elas_2D
        case {4 5}
            SminusA_PEM_elas_2D
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
            ASplus_fluid_fluid_2D
        case 1
            ASplus_elas_fluid_2D
        case {4 5}
            ASplus_PEM_fluid_2D
    end
elseif floor(element_MMT_plus/1000)==1
    switch floor(multilayer(end).mat/1000)
        case {0 2 3}
            fhgfghhghgfgfh
        case 1
            ASplus_elas_elas_2D
        case {4 5}
            ASplus_PEM_elas_2D
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


nb_dof_minus=nS_minus/2;
nb_dof_plus=nS_plus/2;

ZOD_mat=zeros(6,6);
ZOD_mat(index_ZOD,index_ZOD)=T;
% ZOD_mat=T;
