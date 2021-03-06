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



file_PW_id=fopen(name.file_PW,'w');

compute_number_PW_3D



for i_f=1:abs(frequency.nb)
    omega=2*pi*frequency.vec(i_f);
    k_air=omega/air.c;
    k_x=k_air*sin(data_model.theta_inc(1))*cos(data_model.theta_inc(2));
    k_y=k_air*sin(data_model.theta_inc(1))*sin(data_model.theta_inc(2));
    
    if imag(frequency.vec(i_f))==0
    fprintf(file_PW_id,'%1.15e \t',frequency.vec(i_f));
    else
    fprintf(file_PW_id,'%1.15e \t%1.15e \t',real(frequency.vec(i_f)),imag(frequency.vec(i_f)));
    end
    for i_m=1:nb_multilayers
        
        [Mat_PW]=build_global_PW_matrices_3D(k_x,k_y,omega,multilayer(:,i_m),nb_amplitudes(i_m),n_w(:,i_m),k_air,air);
        
        F_PW=-Mat_PW(:,1);
        Mat_PW(:,1)=[];
        
        X_PW=Mat_PW\F_PW;
        
        abs_PW=1-abs(X_PW(1))^2;
        rflx_PW=X_PW(1);
        
        if multilayer(1,1).termination~=0
            trans_PW=X_PW(end);
            TL_PW=-20*log10(abs(X_PW(end)));
            fprintf(file_PW_id,'%1.15e \t%1.15e \t%1.15e \t%1.15e\t%1.15e \t%1.15e \t',abs_PW,real(rflx_PW),imag(rflx_PW),TL_PW,real(trans_PW),imag(trans_PW));
        else
            fprintf(file_PW_id,'%1.15e \t%1.15e \t%1.15e \t%1.15e\t%1.15e \t%1.15e \t',abs_PW,real(rflx_PW),imag(rflx_PW),0,0,0);
        end
        
        
    end
    fprintf(file_PW_id,'\n');
    
end
fclose(file_PW_id);