% PW_resolution.m
%
% Copyright (C) 2016 < Olivier DAZEL <olivier.dazel@univ-lemans.fr> >
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



file_PW_id=fopen(name.file_PW3D,'w');

compute_number_PW_3D



for i_f=1:abs(frequency.nb)
    omega=2*pi*frequency.vec(i_f);
    k_air=omega/air.c;
    k_x=k_air*sin(data_model.theta_1)*cos(data_model.theta_2);
    k_y=k_air*sin(data_model.theta_1)*sin(data_model.theta_2);
    
    fprintf(file_PW_id,'%1.15e \t',frequency.vec(i_f));    
    for i_m=1:nb_multilayers
        
        [Mat_PW,k_z_2,SV_2]=build_global_PW_matrices_3D(k_x,k_y,omega,multilayer(:,i_m),nb_amplitudes(i_m),n_w(:,i_m),k_air,air);

%         medium_2=5010;
%         eval(['Mat_porous_' num2str(medium_2-1000*floor(medium_2/1000))])
%         typ_mat=floor(medium_2/1000);
%         properties_eqf
%         vec_rho(1,i_f)=rho_eq_til(1,1);
%         vec_rho(2,i_f)=rho_eq_til(2,2);
%         vec_rho(3,i_f)=rho_eq_til(3,3);
%         vec_k_z_2(1:4,i_f)=k_z_2;
%         vec_SV_2(i_f)=SV_2(1,3);
%         
        
        
        F_PW=-Mat_PW(:,1);
        Mat_PW(:,1)=[];
        

        X_PW=Mat_PW\F_PW;
        
        abs_PW_3D(i_f,i_m)=1-abs(X_PW(1))^2;
        rflx_PW_3D(i_f,i_m)=X_PW(1);
        vec_X(:,i_f)=X_PW;
        if multilayer(1,1).termination~=0
            TL_PW_3D(i_f,i_m)=-20*log10(abs(X_PW(end,i_m)));
            fprintf(file_PW_id,'%1.15e \t%1.15e \t%1.15e \t%1.15e \t',abs_PW_3D(i_f,i_m),real(rflx_PW_3D(i_f,i_m)),imag(rflx_PW_3D(i_f,i_m)),TL_PW_3D(i_f,i_m));
        else
            fprintf(file_PW_id,'%1.15e \t%1.15e \t%1.15e \t%1.15e \t',abs_PW_3D(i_f,i_m),real(rflx_PW_3D(i_f,i_m)),imag(rflx_PW_3D(i_f,i_m)),0);
        end

        
    end
    fprintf(file_PW_id,'\n');
    
end
fclose(file_PW_id);