

tau_x=exp(j*pi/4);

M_an=[];
F_an=[];


k_0=omega/c_0;

x=-d_pml;
M_an(1,1)=sin(k_0*tau_x*x);
M_an(1,2)=-sin(k_0*x);
M_an(1,3)=-cos(k_0*x);

M_an(2,1)=K_0*k_0*cos(k_0*tau_x*x);
M_an(2,2)=-K_0*k_0*cos(k_0*x);
M_an(2,3)= K_0*k_0*sin(k_0*x);

x=-d_pml-d_air;
M_an(3,2)=sin(k_0*x)*(j*omega);
M_an(3,3)=cos(k_0*x)*(j*omega);

%F_an(2,1)=-1;
F_an(3,1)=-1;

x_poreux=linspace(-(d_pml),0,100);
x_air=linspace(-(d_air+d_pml),-d_pml,100);

X_an=M_an\F_an;

A=X_an(1);
B=X_an(2);
C=X_an(3);

u_x_poreux=A*sin(k_0*tau_x*x_poreux);
v_z_poreux=0*u_x_poreux;
v_x_poreux=j*omega*u_x_poreux;
sigma_p=K_0*k_0*A*cos(k_0*tau_x*x_poreux);

u_x_air=B*sin(k_0*x_air)+C*cos(k_0*x_air);
v_z_air=0*u_x_air;
v_x_air=j*omega*u_x_air;
sigma_a=K_0*k_0*(B*cos(k_0*x_air)-C*sin(k_0*x_air));




% figure(1000)
% subplot(221)
hold on
plot(x_air+d_air,abs(sigma_a))
plot(x_poreux+d_air,abs(sigma_p),'r')
% subplot(223)
% hold on
% plot(x_air,abs(sigma_a))
% plot(x_poreux,abs(sigma_p),'r')
% 
% figure(1001)
% subplot(221)
% hold on
% plot(x_air,angle(v_x_air))
% plot(x_poreux,angle(v_x_poreux),'r')
% subplot(223)
% hold on
% plot(x_air,angle(sigma_a))
% plot(x_poreux,angle(sigma_p),'r')

