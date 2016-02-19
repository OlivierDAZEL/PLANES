% Mat_PEM_3.m
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

porous_model.eqf='jca_aniso';
porous_model.frame='structural';
porous_model.aniso='yes';



phi=0.95;
sig=42000;
alpha=1.100;
LCV=1.50E-05;
LCT=4.500E-05;
rho_1=126.000;
young=694400E+00;
nu=0.24000E+00;
eta=0.05;

N=(young)/(2*(1+nu));
A_hat=(young*nu)/((1+nu)*(1-2*nu));
F_hat=1*A_hat;

C_hat_conservative= [A_hat+2*N A_hat F_hat 0 0 0;A_hat A_hat+2*N F_hat 0 0 0;F_hat F_hat A_hat+2*N 0 0 0; 0 0 0 N 0 0;0 0 0 0 N 0; 0 0 0 0 0 N ];
C_tot=[1 1 1 0 0 0;1 1 1 0 0 0;1 1 1 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0];


%C_hat_conservative=[A_hat+2*N A_hat F_hat 0 0 0;A_hat A_hat+2*N F_hat 0 0 0;F_hat F_hat CC_hat 0 0 0; 0 0 0 L 0 0;0 0 0 0 L 0; 0 0 0 0 0 N ];

% sig_tensor=sig*eye(3);
% alpha_tensor=alpha*eye(3);
% LCV_tensor=LCV*eye(3);

sig_tensor=diag([sig sig sig]);
alpha_tensor=diag([alpha alpha alpha]);
LCV_tensor=diag([LCV LCV LCV]);


alpha_hat=0.33348;
beta_hat =812.69e3;
b_hat=0.29620;

Q=rotate_u3(pi/4,0,0)';

R1 = [ Q(1,1).^2 Q(1,2).^2 Q(1,3).^2 ; ...
       Q(2,1).^2 Q(2,2).^2 Q(2,3).^2 ; ...
       Q(3,1).^2 Q(3,2).^2 Q(3,3).^2 ] ;

R2 = [ Q(1,2).*Q(1,3) Q(1,3).*Q(1,1) Q(1,1).*Q(1,2) ; ...
       Q(2,2).*Q(2,3) Q(2,3).*Q(2,1) Q(2,1).*Q(2,2) ; ...
       Q(3,2).*Q(3,3) Q(3,3).*Q(3,1) Q(3,1).*Q(3,2) ] ;

R3 = [ Q(2,1).*Q(3,1) Q(2,2).*Q(3,2) Q(2,3).*Q(3,3) ; ...
       Q(3,1).*Q(1,1) Q(3,2).*Q(1,2) Q(3,3).*Q(1,3) ; ...
       Q(1,1).*Q(2,1) Q(1,2).*Q(2,2) Q(1,3).*Q(2,3) ] ;

R4 = [ Q(2,2).*Q(3,3)+Q(2,3).*Q(3,2) ...
                     Q(2,3).*Q(3,1)+Q(2,1).*Q(3,3) ...
                                 Q(2,1).*Q(3,2)+Q(2,2).*Q(3,1) ; ...
       Q(3,2).*Q(1,3)+Q(3,3).*Q(1,2) ...
                     Q(3,3).*Q(1,1)+Q(3,1).*Q(1,3) ...
                                 Q(3,1).*Q(1,2)+Q(3,2).*Q(1,1) ; ...      
       Q(1,2).*Q(2,3)+Q(1,3).*Q(2,2) ...
                     Q(1,3).*Q(2,1)+Q(1,1).*Q(2,3) ...
                                 Q(1,1).*Q(2,2)+Q(1,2).*Q(2,1)] ; 

M_rot = [ R1  2*R2 ; R3   R4   ] ;

C_hat_conservative=M_rot * C_hat_conservative * transpose(M_rot) ;
C_tot=M_rot * C_tot * transpose(M_rot) ;



% % % % % % % % Mat_PEM_3.m
% % % % % % % %
% % % % % % % % Copyright (C) 2014 < Olivier DAZEL <olivier.dazel@univ-lemans.fr> >
% % % % % % % %
% % % % % % % % This file is part of PLANES.
% % % % % % % %
% % % % % % % % PLANES (Porous LAum NumErical Simulator) is a software to compute the
% % % % % % % % vibroacoustic response of sound packages containing coupled
% % % % % % % % acoustic/elastic/porous substructures. It is mainly based on the
% % % % % % % % Finite-Element Method and some numerical methods developped at
% % % % % % % % LAUM (http://laum.univ-lemans.fr).
% % % % % % % %
% % % % % % % % You can download the latest version of PLANES at
% % % % % % % % https://github.com/OlivierDAZEL/PLANES
% % % % % % % % or find more details on Olivier's webpage
% % % % % % % % http://perso.univ-lemans.fr/~odazel/
% % % % % % % %
% % % % % % % % For any questions or if you want to
% % % % % % % % contribute to this project, contact Olivier.
% % % % % % % %
% % % % % % % % PLANES is free software: you can redistribute it and/or modify
% % % % % % % % it under the terms of the GNU General Public License as published by
% % % % % % % % the Free Software Foundation, either version 3 of the License, or
% % % % % % % % (at your option) any later version.
% % % % % % % %
% % % % % % % % This program is distributed in the hope that it will be useful,
% % % % % % % % but WITHOUT ANY WARRANTY; without even the implied warranty of
% % % % % % % % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% % % % % % % % GNU General Public License for more details.
% % % % % % % %
% % % % % % % % You should have received a copy of the GNU General Public License
% % % % % % % % along with this program. If not, see <http://www.gnu.org/licenses/>.
% % % % % % % %%
% % % % % % % 
% % % % % % % porous_model.eqf='jca_aniso';
% % % % % % % porous_model.frame='structural';
% % % % % % % porous_model.aniso='yes';
% % % % % % % 
% % % % % % % 
% % % % % % % % phi=0.99;
% % % % % % % %
% % % % % % % % sig_z=8000;
% % % % % % % % sig_x=4000;
% % % % % % % %
% % % % % % % % alpha_z=1.1;
% % % % % % % % alpha_x=1.1;
% % % % % % % %
% % % % % % % % LCV_z=140*1e-6;
% % % % % % % % LCV_x=200*1e-6;
% % % % % % % %
% % % % % % % % LCT=500*1e-6;
% % % % % % % %
% % % % % % % % rho_1=32;
% % % % % % % % eta=0.01;
% % % % % % % %
% % % % % % % %
% % % % % % % % G=260*1e3;
% % % % % % % % Gp=125*1e3;
% % % % % % % % CC_hat=46*1e3;
% % % % % % % %
% % % % % % % % N=G;
% % % % % % % % L=Gp;
% % % % % % % % A_hat=0;
% % % % % % % % F_hat=0;
% % % % % % % 
% % % % % % % 
% % % % % % % C_hat_conservative=[A_hat+2*N A_hat F_hat 0 0 0;A_hat A_hat+2*N F_hat 0 0 0;F_hat F_hat CC_hat 0 0 0; 0 0 0 L 0 0;0 0 0 0 L 0; 0 0 0 0 0 N ];
% % % % % % % 
% % % % % % % % sig_tensor=sig*eye(3);
% % % % % % % % alpha_tensor=alpha*eye(3);
% % % % % % % % LCV_tensor=LCV*eye(3);
% % % % % % % 
% % % % % % % sig_tensor=diag([sig_x sig_x sig_z]);
% % % % % % % alpha_tensor=diag([alpha_x alpha_x alpha_z]);
% % % % % % % LCV_tensor=diag([LCV_x LCV_x LCV_z]);
% % % % % % % 
% % % % % % % 
% % % % % % % alpha_hat=0.33348;
% % % % % % % beta_hat =812.69e3;
% % % % % % % b_hat=0.29620;
% % % % % % % 
% % % % % % % Q=rotate_u3(0,0,0)';
% % % % % % % 
% % % % % % % R1 = [ Q(1,1).^2 Q(1,2).^2 Q(1,3).^2 ; ...
% % % % % % %     Q(2,1).^2 Q(2,2).^2 Q(2,3).^2 ; ...
% % % % % % %     Q(3,1).^2 Q(3,2).^2 Q(3,3).^2 ] ;
% % % % % % % 
% % % % % % % R2 = [ Q(1,2).*Q(1,3) Q(1,3).*Q(1,1) Q(1,1).*Q(1,2) ; ...
% % % % % % %     Q(2,2).*Q(2,3) Q(2,3).*Q(2,1) Q(2,1).*Q(2,2) ; ...
% % % % % % %     Q(3,2).*Q(3,3) Q(3,3).*Q(3,1) Q(3,1).*Q(3,2) ] ;
% % % % % % % 
% % % % % % % R3 = [ Q(2,1).*Q(3,1) Q(2,2).*Q(3,2) Q(2,3).*Q(3,3) ; ...
% % % % % % %     Q(3,1).*Q(1,1) Q(3,2).*Q(1,2) Q(3,3).*Q(1,3) ; ...
% % % % % % %     Q(1,1).*Q(2,1) Q(1,2).*Q(2,2) Q(1,3).*Q(2,3) ] ;
% % % % % % % 
% % % % % % % R4 = [ Q(2,2).*Q(3,3)+Q(2,3).*Q(3,2) ...
% % % % % % %     Q(2,3).*Q(3,1)+Q(2,1).*Q(3,3) ...
% % % % % % %     Q(2,1).*Q(3,2)+Q(2,2).*Q(3,1) ; ...
% % % % % % %     Q(3,2).*Q(1,3)+Q(3,3).*Q(1,2) ...
% % % % % % %     Q(3,3).*Q(1,1)+Q(3,1).*Q(1,3) ...
% % % % % % %     Q(3,1).*Q(1,2)+Q(3,2).*Q(1,1) ; ...
% % % % % % %     Q(1,2).*Q(2,3)+Q(1,3).*Q(2,2) ...
% % % % % % %     Q(1,3).*Q(2,1)+Q(1,1).*Q(2,3) ...
% % % % % % %     Q(1,1).*Q(2,2)+Q(1,2).*Q(2,1)] ;
% % % % % % % 
% % % % % % % M_rot = [ R1  2*R2 ; R3   R4   ] ;
% % % % % % % 
% % % % % % % C_hat_conservative=M_rot * C_hat_conservative * transpose(M_rot) ;