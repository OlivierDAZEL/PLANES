% solution_air_EF.m
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

% Dans le poreux
%% u^t= mu_1 *A sin(delta_1*x)+mu_2*B sin(delta_2*x)
%% sigma_p=K_eq_til*(mu_1*delta_1*A*cos(delta_1*x)+mu_2*delta_2*B*cos(delta_2*x));
% Dans l'air
%% u^a= C sin(k_0*x)+D cos(k_0*x)
%% sigma_a=K_0*k_0(C*cos(k_0*x)-D*sin(k_0*x)

if solve.DGM
rho_eq_til=Mat_parameter(1,2);
K_eq_til=Mat_parameter(2,2);
end

delta_eq=omega*sqrt(rho_eq_til/K_eq_til);

M_an=zeros(3,3);
F_an=zeros(3,1);

x=-d_EF;

% Continuite du d?placement
M_an(1,1)=sin(delta_eq*x);
M_an(1,2)=-sin(k_air*x);
M_an(1,3)=-cos(k_air*x);
% Continuite de la pression
M_an(2,1)=K_eq_til*delta_eq*cos(delta_eq*x);
M_an(2,2)=-air.K*k_air*cos(k_air*x);
M_an(2,3)= air.K*k_air*sin(k_air*x);

x=-d_EF-d_air;
% % Pression unite

% Unite displacement

M_an(3,2)= sin(k_air*x);
M_an(3,3)= cos(k_air*x);
F_an(3,1)=-1;




X_an=M_an\F_an;

A=X_an(1);
C=X_an(2);
D=X_an(3);

x_EF=linspace(-(d_EF),0,1000);
x_air=linspace(-(d_air+d_EF),-d_EF,1000);


sigma_p=K_eq_til*delta_eq*A*cos(delta_eq*x_EF);
sigma_a=air.K*k_air*(C*cos(k_air*x_air)-D*sin(k_air*x_air));
%sigma_a=K_eq_til*delta_eq*(C*cos(delta_eq*x_air)-D*sin(delta_eq*x_air));


if plot_profiles==1
    
    
    % figure (1001)
    % hold on
    
    
    
    figure(10002)
    hold on
    plot(x_air+d_air+d_EF,abs(sigma_a))
    plot(x_EF+d_air+d_EF,abs(sigma_p),'m')
    
    % figure(1000)
    % subplot(3,3,1)
    % plot(x_EF,abs(v_s),'r')
    % subplot(3,3,3)
    % hold on
    % plot(x_EF,abs(v_t),'r')
    % plot(x_air,abs(v_x_air))
    % subplot(3,3,5)
    % plot(x_EF,abs(sigma_xx+sigma_yy)/2,'r')
    % subplot(3,3,8)
    % hold on
    % plot(x_EF,abs(sigma_p),'r')
    % plot(x_air,abs(sigma_a))
    % subplot(3,3,2)
    % plot(x_EF,abs(0),'r')
    % subplot(3,3,4)
    % hold on
    % plot(x_air,abs(0))
    % plot(x_EF,abs(0),'r')
    % subplot(3,3,6)
    % plot(x_EF,abs(0),'r')
    % subplot(3,3,7)
    % plot(x_EF,abs(sigma_xx-sigma_yy)/2,'r')
    %
    % %
    % %
    % %
    % % % figure(1001)
    % % % subplot(3,3,1)
    % % % plot(x_EF+d_poreux+d_air,angle(v_s),'r')
    % % % subplot(3,3,3)
    % % % hold on
    % % % plot(x_EF+d_poreux+d_air,angle(v_t),'r')
    % % % plot(x_air+d_poreux+d_air,angle(v_x_air))
    % % % subplot(3,3,5)
    % % % plot(x_EF+d_poreux+d_air,angle(sigma_xx+sigma_yy),'r')
    % % % subplot(3,3,8)
    % % % hold on
    % % % plot(x_EF+d_poreux+d_air,angle(sigma_p),'r')
    % % % plot(x_air+d_poreux+d_air,angle(sigma_a))
    % % % subplot(3,3,2)
    % % % plot(x_EF+d_poreux+d_air,angle(0),'r')
    % % % subplot(3,3,4)
    % % % hold on
    % % % plot(x_air+d_poreux+d_air,angle(0))
    % % % plot(x_EF+d_poreux+d_air,angle(0),'r')
    % % % subplot(3,3,6)
    % % % plot(x_EF+d_poreux+d_air,angle(0),'r')
    % % % subplot(3,3,7)
    % % % plot(x_EF+d_poreux+d_air,angle(sigma_xx-sigma_yy),'r')
    % %
end



% for e_edge=1:nb_elements
%     centre_elem=mean(nodes(elements(e_edge,:),1:2))';
%     parameter_element
%     e_plus=exp(+1j*k_e*(centre_elem(2)-(d_air+d_EF)));
%     e_moins=1/e_plus;
%     
%     if floor(element_label(e_edge)/1000)==2
%           q_theorique(dof_start_element(e_edge)+1)=e_plus*omega*A/2;
%         q_theorique(dof_start_element(e_edge)  )=e_moins*omega*A/2;
%     end
%     if  floor(element_label(e_edge)/1000)==0
%         
% 
%         q_theorique(dof_start_element(e_edge)+1)=e_plus*omega*((C/2-D/(2*i)));
%         q_theorique(dof_start_element(e_edge)  )=e_moins*omega*((C/2+D/(2*i)));
%         
%     end
% end
% 
% 
% figure (21)
% hold on
% 
% 
% plot(abs(X),'r.')
% plot(abs(q_theorique),'k+')
%plot(abs(q_theorique*X(27)/q_theorique(27)),'k+')
% 
% figure (22)
% hold on
% plot(angle(X),'r.')
% plot(angle(-q_theorique),'k+')
% 
