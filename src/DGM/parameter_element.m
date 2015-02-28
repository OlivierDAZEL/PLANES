% parameter_element.m
%
% Copyright (C) 2015 < Olivier DAZEL <olivier.dazel@univ-lemans.fr> >
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
%                 https://github.com/OlivierDAZEL/PLANES
% or find more details on Olivier's webpage
%                  http://perso.univ-lemans.fr/~odazel/
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
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

if floor(element_label(e_edge)/1000)==0
    c_e=air.c;
    k_e=omega/c_e;
    Z_e=air.Z;
    M_e=diag([air.rho,air.rho,1/air.K]);

    tau_x=1;
    tau_y=1;
    
    
elseif floor(element_label(e_edge)/1000)==2
    
    tau_x=1;
    tau_y=1;
        
    k_e=omega*sqrt(Mat_parameter(1,index_element(e_edge))/Mat_parameter(2,index_element(e_edge)));
    c_e=omega/k_e;
    Z_e=Mat_parameter(1,index_element(e_edge))*c_e;
    M_e=diag([Mat_parameter(1,index_element(e_edge)),Mat_parameter(1,index_element(e_edge)),1/Mat_parameter(2,index_element(e_edge))]); 
    

elseif floor(element_label(e_edge)/1000)==8    
    c_e=air.c;
    
    tau_x=1;
    tau_y=1;
    temp=(element_label(e_edge)-8000);
    if (floor(temp/100)==1)
        tau_x=exp(1j*pi/4);
    end
    temp=temp-100*floor(temp/100);
    if (floor(temp/10)==1)
        tau_y=exp(1j*pi/4);
    end    
    
    k_e=omega/c_e;
    Z_e=air.Z;
    
    M_e=diag([air.rho,air.rho,1/air.K]);    
    
elseif floor(element_label(e_edge)/1000)==4
    
    tau_x=1;
    tau_y=1;
    
    rho_eq_til=Mat_parameter(1,index_element(e_edge));
    K_eq_til  =Mat_parameter(2,index_element(e_edge));
    gamma_til =Mat_parameter(3,index_element(e_edge));
    A_hat     =Mat_parameter(4,index_element(e_edge));
    P_hat     =Mat_parameter(5,index_element(e_edge));
    N         =Mat_parameter(6,index_element(e_edge));
    rho_s_til =Mat_parameter(7,index_element(e_edge));
    rho_til   =Mat_parameter(8,index_element(e_edge));
    phi       =Mat_parameter(9,index_element(e_edge));
    compute_Biot_waves

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