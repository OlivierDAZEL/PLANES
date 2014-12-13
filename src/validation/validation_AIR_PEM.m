

d_y=d_PEM;
calcul_propriete
calcul_theorique

k_0=omega/c_0;

% Dans le poreux
%% u^s= A sin(delta_1*x)+B sin(delta_2*x)
%% u^t= mu_1 *A sin(delta_1*x)+mu_2*B sin(delta_2*x)
%% sigma_xx=P_hat*(delta_1*A*cos(delta_1*x)+delta_2*B*cos(delta_2*x));
%% sigma_p=K_eq_til*(mu_1*delta_1*A*cos(delta_1*x)+mu_2*delta_2*B*cos(delta_2*x));
% Dans l'air
%% u^a= C sin(k_0*x)+D cos(k_0*x)
%% sigma_a=K_0*k_0(C*cos(k_0*x)-D*sin(k_0*x)


M_an=zeros(4,4);
F_an=zeros(4,1);

x=-d_PEM;

% Continuite du d?placement
M_an(1,1)=mu_1*sin(delta_1*x);
M_an(1,2)=mu_2*sin(delta_2*x);
M_an(1,3)=-sin(k_0*x);
M_an(1,4)=-cos(k_0*x);
% Continuite de la pression
M_an(2,1)=K_eq_til*delta_1*mu_1*cos(delta_1*x);
M_an(2,2)=K_eq_til*delta_2*mu_2*cos(delta_2*x);
M_an(2,3)=-K_0*k_0*cos(k_0*x);
M_an(2,4)= K_0*k_0*sin(k_0*x);
% Nullit? in vacuo
M_an(3,1)=delta_1*cos(delta_1*x);
M_an(3,2)=delta_2*cos(delta_2*x);

x=-d_PEM-d_air;
% % Pression unite
%
% M_an(4,3)= K_0*k_0*cos(k_0*x);
% M_an(4,4)=-K_0*k_0*sin(k_0*x);
% F_an(4,1)=-1;


% Velocity unite

M_an(4,3)= sin(k_0*x);
M_an(4,4)= cos(k_0*x);
F_an(4,1)=-1/(j*omega);




X_an=M_an\F_an;

A=X_an(1);
B=X_an(2);
C=X_an(3);
D=X_an(4);

x_poreux=linspace(-(d_PEM),0,100);
x_air=linspace(-(d_air+d_PEM),-d_PEM,100);

u_s=A*sin(delta_1*x_poreux)+B*sin(delta_2*x_poreux);
u_t=A*mu_1*sin(delta_1*x_poreux)+B*mu_2*sin(delta_2*x_poreux);
v_s=j*omega*u_s;
v_t=j*omega*u_t;
sigma_xx=P_hat*(delta_1*A*cos(delta_1*x_poreux)+delta_2*B*cos(delta_2*x_poreux));
sigma_yy=A_hat*(delta_1*A*cos(delta_1*x_poreux)+delta_2*B*cos(delta_2*x_poreux));
sigma_p=K_eq_til*(mu_1*delta_1*A*cos(delta_1*x_poreux)+mu_2*delta_2*B*cos(delta_2*x_poreux));


u_x_air=C*sin(k_0*x_air)+D*cos(k_0*x_air);
v_z_air=0*u_x_air;
v_x_air=j*omega*u_x_air;
sigma_a=K_0*k_0*(C*cos(k_0*x_air)-D*sin(k_0*x_air));


% % Calcul des normes L2
% 
% 
% L2_norm_analytique_ua=conj(C)*C*SS(k_0,k_0,-(d_air+d_PEM),-d_PEM)+...
%                       conj(C)*D*SC(k_0,k_0,-(d_air+d_PEM),-d_PEM)+...
%                       conj(D)*C*CS(k_0,k_0,-(d_air+d_PEM),-d_PEM)+...
%                       conj(D)*D*CC(k_0,k_0,-(d_air+d_PEM),-d_PEM);
%                  L2_norm_analytique_ua=sqrt(real(L2_norm_analytique_ua))*sqrt(d_y);
%                  
%                  
% % valid_L2_ua=norm(u_x_air)^2-(norm(u_x_air(1))^2+norm(u_x_air(end))^2)/2;
% % valid_L2_ua=valid_L2_ua*(d_air)/(length(u_x_air)-1);
% % valid_L2_ua=sqrt(real(valid_L2_ua));
%                  
%           
% L2_norm_analytique_sigma_a= conj(C)*C*CC(k_0,k_0,-(d_air+d_PEM),-d_PEM)...
%                            -conj(C)*D*CS(k_0,k_0,-(d_air+d_PEM),-d_PEM)...
%                            -conj(D)*C*SC(k_0,k_0,-(d_air+d_PEM),-d_PEM)+...
%                             conj(D)*D*SS(k_0,k_0,-(d_air+d_PEM),-d_PEM);
% L2_norm_analytique_sigma_a=omega*norm(K_0*k_0)*sqrt(L2_norm_analytique_sigma_a)*sqrt(d_y);
% 
% 
% % valid_L2_sigma_a=norm(sigma_a)^2-(norm(sigma_a(1))^2+norm(sigma_a(end))^2)/2;
% % valid_L2_sigma_a=valid_L2_sigma_a*(d_air)/(length(u_x_air)-1);
% % valid_L2_sigma_a=sqrt(real(valid_L2_sigma_a));
% 
% 
% L2_norm_analytique_us=conj(A)*A*SS(delta_1,delta_1,-d_PEM,0)+...
%                       conj(A)*B*SS(delta_1,delta_2,-d_PEM,0)+...
%                       conj(B)*A*SS(delta_2,delta_1,-d_PEM,0)+...
%                       conj(B)*B*SS(delta_2,delta_2,-d_PEM,0);
% L2_norm_analytique_us=omega*sqrt(real(L2_norm_analytique_us))*sqrt(d_y);
%                  
%                  
% % valid_L2_us=norm(u_s)^2-(norm(u_s(1))^2+norm(u_s(end))^2)/2;
% % valid_L2_us=valid_L2_us*(d_PEM)/(length(u_s)-1);
% % valid_L2_us=sqrt(real(valid_L2_us))*d_y
%                  
%           
% L2_norm_analytique_sigma_p=conj(mu_1*delta_1*A)*mu_1*delta_1*A*CC(delta_1,delta_1,-d_PEM,0)+...
%                            conj(mu_1*delta_1*A)*mu_2*delta_2*B*CC(delta_1,delta_2,-d_PEM,0)+...
%                            conj(mu_2*delta_2*B)*mu_1*delta_1*A*CC(delta_2,delta_1,-d_PEM,0)+...
%                            conj(mu_2*delta_2*B)*mu_2*delta_2*B*CC(delta_2,delta_2,-d_PEM,0);
% 
% L2_norm_analytique_sigma_p=omega*norm(K_eq_til)*sqrt(L2_norm_analytique_sigma_p)*sqrt(d_y);
%                   
%                  
%                   
%                   
% 
% % valid_L2_sigma_p=norm(sigma_p)^2-(norm(sigma_p(1))^2+norm(sigma_p(end))^2)/2;
% % valid_L2_sigma_p=valid_L2_sigma_p*(d_PEM)/(length(u_s)-1);
% % valid_L2_sigma_p=sqrt(real(valid_L2_sigma_p))*d_y
% 
%           
% 
% % if tracefigure==1
% % 
% %     figure (1000)
% %     set(gca,'Fontsize',15)
% %     hold on
% % 
% %     figure(1000)
% %     plot(x_poreux,abs(v_s),'r','Linewidth',4)
% %     figure (1001)
% %     set(gca,'Fontsize',15)
% %     hold on
% % 
% %     plot(x_poreux,abs(sigma_p),'r','Linewidth',4)
% %     plot(x_air,abs(sigma_a),'Linewidth',4)
% %     
% % 
% % end














% figure (1000)
% hold on
% % figure (1001)
% % hold on
% 
% figure(1000)
% subplot(3,3,1)
% plot(x_poreux,abs(v_s),'r')
% subplot(3,3,3)
hold on
plot(x_poreux+d_air,abs(v_t),'m')
plot(x_air+d_air,abs(v_x_air),'m')
% subplot(3,3,5)
% plot(x_poreux,abs(sigma_xx+sigma_yy)/2,'r')
% subplot(3,3,8)
%hold on
%plot(x_poreux,abs(sigma_p),'r')
%plot(x_air,abs(sigma_a))
%subplot(3,3,2)
% plot(x_poreux,abs(0),'r')
% subplot(3,3,4)
% hold on
% plot(x_air,abs(0))
% plot(x_poreux,abs(0),'r')
% subplot(3,3,6)
% plot(x_poreux,abs(0),'r')
% subplot(3,3,7)
% plot(x_poreux,abs(sigma_xx-sigma_yy)/2,'r')
% 
% %
% %
% %
% % % figure(1001)
% % % subplot(3,3,1)
% % % plot(x_poreux+d_PEM+d_air,angle(v_s),'r')
% % % subplot(3,3,3)
% % % hold on
% % % plot(x_poreux+d_PEM+d_air,angle(v_t),'r')
% % % plot(x_air+d_PEM+d_air,angle(v_x_air))
% % % subplot(3,3,5)
% % % plot(x_poreux+d_PEM+d_air,angle(sigma_xx+sigma_yy),'r')
% % % subplot(3,3,8)
% % % hold on
% % % plot(x_poreux+d_PEM+d_air,angle(sigma_p),'r')
% % % plot(x_air+d_PEM+d_air,angle(sigma_a))
% % % subplot(3,3,2)
% % % plot(x_poreux+d_PEM+d_air,angle(0),'r')
% % % subplot(3,3,4)
% % % hold on
% % % plot(x_air+d_PEM+d_air,angle(0))
% % % plot(x_poreux+d_PEM+d_air,angle(0),'r')
% % % subplot(3,3,6)
% % % plot(x_poreux+d_PEM+d_air,angle(0),'r')
% % % subplot(3,3,7)
% % % plot(x_poreux+d_PEM+d_air,angle(sigma_xx-sigma_yy),'r')
% %
% 
% 
% 
% 
% 
% 
