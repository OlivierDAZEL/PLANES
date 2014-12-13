function M=State_PEM(k_x,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til)



beta_1=sqrt(delta_1^2-k_x^2);
beta_2=sqrt(delta_2^2-k_x^2);
beta_3=sqrt(delta_3^2-k_x^2);

alpha_1=-1i*A_hat*delta_1^2-1i*2*N*beta_1^2;
alpha_2=-1i*A_hat*delta_2^2-1i*2*N*beta_2^2;
alpha_3= 2*1i*N*beta_3*k_x;


M(1:6,1)=[-2*1i*N*beta_1*k_x; beta_1; mu_1*beta_1;alpha_1;1i*delta_1^2*K_eq_til*mu_1;k_x];
M(1:6,4)=[ 2*1i*N*beta_1*k_x;-beta_1;-mu_1*beta_1;alpha_1;1i*delta_1^2*K_eq_til*mu_1;k_x];

M(1:6,2)=[-2*1i*N*beta_2*k_x; beta_2; mu_2*beta_2;alpha_2;1i*delta_2^2*K_eq_til*mu_2;k_x];
M(1:6,5)=[ 2*1i*N*beta_2*k_x;-beta_2;-mu_2*beta_2;alpha_2;1i*delta_2^2*K_eq_til*mu_2;k_x];

M(1:6,3)=[1i*N*(beta_3^2-k_x^2);k_x;mu_3*k_x;-alpha_3;0;-beta_3];
M(1:6,6)=[1i*N*(beta_3^2-k_x^2);k_x;mu_3*k_x; alpha_3;0; beta_3];