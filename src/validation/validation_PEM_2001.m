


calcul_propriete;


calcul_theorique;

k_0=omega/c_0;
d_poreux=dyair;


% Dans le poreux
%% u^s= A sin(delta_1*x)+B sin(delta_2*x)
%% u^t= mu_1 *A sin(delta_1*x)+mu_2*B sin(delta_2*x)
%% sigma_xx=P_hat*(delta_1*A*cos(delta_1*x)+delta_2*B*cos(delta_2*x));
%% sigma_p=K_eq_til*(mu_1*delta_1*A*cos(delta_1*x)+mu_2*delta_2*B*cos(delta_2*x));


M_an=zeros(2,2);
F_an=zeros(2,1);



x=-d_poreux;
% sigma_total=1
M_an(1,1)=P_hat*(delta_1*cos(delta_1*x))+K_eq_til*delta_1*mu_1*cos(delta_1*x);
M_an(1,2)=P_hat*(delta_2*cos(delta_2*x))+K_eq_til*delta_2*mu_2*cos(delta_2*x);
F_an(1,1)=1;

% Nullit? deplacment relatif
M_an(2,1)=(mu_1-1)*sin(delta_1*x);
M_an(2,2)=(mu_2-1)*sin(delta_2*x);



X_an=M_an\F_an;

A=X_an(1);
B=X_an(2);

x_poreux=linspace(-(d_poreux),0,100);


u_s=A*sin(delta_1*x_poreux)+B*sin(delta_2*x_poreux);
u_t=A*mu_1*sin(delta_1*x_poreux)+B*mu_2*sin(delta_2*x_poreux);
v_s=j*omega*u_s;
v_t=j*omega*u_t;
sigma_xx=P_hat*(delta_1*A*cos(delta_1*x_poreux)+delta_2*B*cos(delta_2*x_poreux));
sigma_yy=A_hat*(delta_1*A*cos(delta_1*x_poreux)+delta_2*B*cos(delta_2*x_poreux));
sigma_p=K_eq_til*(mu_1*delta_1*A*cos(delta_1*x_poreux)+mu_2*delta_2*B*cos(delta_2*x_poreux));


d_x_air=0;
d_x_poreux=d_poreux;



figure(5678)
hold on
plot(x_poreux+d_x_poreux+d_x_air,abs(u_s),'r')
plot(x_poreux+d_x_poreux+d_x_air,abs(u_t),'m')

figure(5679)
plot(x_poreux+d_x_poreux+d_x_air,abs(-sigma_p),'r')

 Z_analytique=(-sigma_p(1))/(v_t(1));
 R_analytique=(Z_analytique-Z_0)/(Z_analytique+Z_0);

abs_analytique(i_f)=1-abs(R_analytique)^2;



% figure(56780)
% plot(x_poreux+d_x_poreux+d_x_air,angle(u_s),'r')
% figure(56790)
% plot(x_poreux+d_x_poreux+d_x_air,angle(-sigma_p),'r')
% subplot(3,3,3)
% hold on
% plot(x_poreux+d_x_poreux+d_x_air,abs(v_t),'r')
% subplot(3,3,5)
% plot(x_poreux+d_x_poreux+d_x_air,abs(sigma_xx+sigma_yy)/2,'r')
% subplot(3,3,8)
% hold on
% plot(x_poreux+d_x_poreux+d_x_air,abs(sigma_p),'r')
% subplot(3,3,2)
% plot(x_poreux+d_x_poreux+d_x_air,abs(0),'r')
% subplot(3,3,4)
% hold on
% plot(x_poreux+d_x_poreux+d_x_air,abs(0),'r')
% subplot(3,3,6)
% plot(x_poreux+d_x_poreux+d_x_air,abs(0),'r')
% subplot(3,3,7)
% plot(x_poreux+d_x_poreux+d_x_air,abs(sigma_xx-sigma_yy)/2,'r')
% 
% 
% 
% 
% figure(1001)
% subplot(3,3,1)
% plot(x_poreux+d_x_poreux+d_x_air,angle(v_s),'r')
% subplot(3,3,3)
% hold on
% plot(x_poreux+d_x_poreux+d_x_air,angle(v_t),'r')
% subplot(3,3,5)
% plot(x_poreux+d_x_poreux+d_x_air,angle(sigma_xx+sigma_yy),'r')
% subplot(3,3,8)
% hold on
% plot(x_poreux+d_x_poreux+d_x_air,angle(sigma_p),'r')
% subplot(3,3,2)
% plot(x_poreux+d_x_poreux+d_x_air,angle(0),'r')
% subplot(3,3,4)
% hold on
% plot(x_poreux+d_x_poreux+d_x_air,angle(0),'r')
% subplot(3,3,6)
% plot(x_poreux+d_x_poreux+d_x_air,angle(0),'r')
% subplot(3,3,7)
% plot(x_poreux+d_x_poreux+d_x_air,angle(sigma_xx-sigma_yy),'r')







