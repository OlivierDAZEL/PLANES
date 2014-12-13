function M=State_elas(k_x,k_y,delta_P,delta_s,lambda,mu)


beta_P=sqrt(delta_P^2-k_x^2);
beta_s=sqrt(delta_s^2-k_x^2);


alpha_P=-j*lambda*delta_P^2-j*2*mu*beta_P^2;
alpha_s= 2*j*mu*beta_s*k_x;


M(1:4,1)=[-2*j*mu*beta_P*k_x; beta_P;alpha_P;k_x];
M(1:4,3)=[ 2*j*mu*beta_P*k_x;-beta_P;alpha_P;k_x];

M(1:4,2)=[j*mu*(beta_s^2-k_x^2);k_x;-alpha_s;-beta_s];
M(1:4,4)=[j*mu*(beta_s^2-k_x^2);k_x; alpha_s; beta_s];
