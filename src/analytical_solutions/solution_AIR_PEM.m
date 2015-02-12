
compute_Biot_waves

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
M_an(1,3)=-sin(k_air*x);
M_an(1,4)=-cos(k_air*x);
% Continuite de la pression
M_an(2,1)=K_eq_til*delta_1*mu_1*cos(delta_1*x);
M_an(2,2)=K_eq_til*delta_2*mu_2*cos(delta_2*x);
M_an(2,3)=-air.K*k_air*cos(k_air*x);
M_an(2,4)= air.K*k_air*sin(k_air*x);
% Nullit? in vacuo
M_an(3,1)=delta_1*cos(delta_1*x);
M_an(3,2)=delta_2*cos(delta_2*x);

x=-d_PEM-d_air;
% % Pression unite
%
% M_an(4,3)= K_0*k_0*cos(k_0*x);
% M_an(4,4)=-K_0*k_0*sin(k_0*x);
% F_an(4,1)=-1;


% Unite displacement

M_an(4,3)= sin(k_air*x);
M_an(4,4)= cos(k_air*x);
F_an(4,1)=-1;




X_an=M_an\F_an;

A=X_an(1);
B=X_an(2);
C=X_an(3);
D=X_an(4);

x_PEM=linspace(-(d_PEM),0,100);
x_air=linspace(-(d_air+d_PEM),-d_PEM,100);

u_s=A*sin(delta_1*x_PEM)+B*sin(delta_2*x_PEM);
u_t=A*mu_1*sin(delta_1*x_PEM)+B*mu_2*sin(delta_2*x_PEM);
v_s=j*omega*u_s;
v_t=j*omega*u_t;
sigma_xx=P_hat*(delta_1*A*cos(delta_1*x_PEM)+delta_2*B*cos(delta_2*x_PEM));
sigma_yy=A_hat*(delta_1*A*cos(delta_1*x_PEM)+delta_2*B*cos(delta_2*x_PEM));
sigma_p=K_eq_til*(mu_1*delta_1*A*cos(delta_1*x_PEM)+mu_2*delta_2*B*cos(delta_2*x_PEM));


u_x_air=C*sin(k_air*x_air)+D*cos(k_air*x_air);
v_z_air=0*u_x_air;
v_x_air=j*omega*u_x_air;
sigma_a=air.K*k_air*(C*cos(k_air*x_air)-D*sin(k_air*x_air));




if plot_profiles==1
    

    % figure (1001)
    % hold on
    
    
    figure(10001)
    plot(x_PEM+d_air+d_PEM,abs(u_s),'k')
    
    figure(10002)
    hold on
    plot(x_air+d_air+d_PEM,abs(sigma_a))
    plot(x_PEM+d_air+d_PEM,abs(sigma_p),'m')
    
    % figure(1000)
    % subplot(3,3,1)
    % plot(x_PEM,abs(v_s),'r')
    % subplot(3,3,3)
    % hold on
    % plot(x_PEM,abs(v_t),'r')
    % plot(x_air,abs(v_x_air))
    % subplot(3,3,5)
    % plot(x_PEM,abs(sigma_xx+sigma_yy)/2,'r')
    % subplot(3,3,8)
    % hold on
    % plot(x_PEM,abs(sigma_p),'r')
    % plot(x_air,abs(sigma_a))
    % subplot(3,3,2)
    % plot(x_PEM,abs(0),'r')
    % subplot(3,3,4)
    % hold on
    % plot(x_air,abs(0))
    % plot(x_PEM,abs(0),'r')
    % subplot(3,3,6)
    % plot(x_PEM,abs(0),'r')
    % subplot(3,3,7)
    % plot(x_PEM,abs(sigma_xx-sigma_yy)/2,'r')
    %
    % %
    % %
    % %
    % % % figure(1001)
    % % % subplot(3,3,1)
    % % % plot(x_PEM+d_poreux+d_air,angle(v_s),'r')
    % % % subplot(3,3,3)
    % % % hold on
    % % % plot(x_PEM+d_poreux+d_air,angle(v_t),'r')
    % % % plot(x_air+d_poreux+d_air,angle(v_x_air))
    % % % subplot(3,3,5)
    % % % plot(x_PEM+d_poreux+d_air,angle(sigma_xx+sigma_yy),'r')
    % % % subplot(3,3,8)
    % % % hold on
    % % % plot(x_PEM+d_poreux+d_air,angle(sigma_p),'r')
    % % % plot(x_air+d_poreux+d_air,angle(sigma_a))
    % % % subplot(3,3,2)
    % % % plot(x_PEM+d_poreux+d_air,angle(0),'r')
    % % % subplot(3,3,4)
    % % % hold on
    % % % plot(x_air+d_poreux+d_air,angle(0))
    % % % plot(x_PEM+d_poreux+d_air,angle(0),'r')
    % % % subplot(3,3,6)
    % % % plot(x_PEM+d_poreux+d_air,angle(0),'r')
    % % % subplot(3,3,7)
    % % % plot(x_PEM+d_poreux+d_air,angle(sigma_xx-sigma_yy),'r')
    % %
end






