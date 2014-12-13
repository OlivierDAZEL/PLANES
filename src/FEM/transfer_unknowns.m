function    [Omega_plus,T_back]=transfer_unknowns(k_x,omega,Omega_moins,signe,media)



lambda=media(1).lambda;
mu=media(1).mu;
rho=media(1).rho;
d=signe*media(1).thickness;

P=lambda+2*mu;


delta_P=omega*sqrt(rho/(P));
delta_s=omega*sqrt(rho/(mu));


beta_P=sqrt(delta_P^2-k_x^2);
beta_s=sqrt(delta_s^2-k_x^2);


alpha_P=-1j*lambda*delta_P^2-1j*2*mu*beta_P^2;
alpha_s= 2*1i*mu*beta_s*k_x;

 V_0=diag([1i*beta_P,-1i*beta_P,1i*beta_s,-1i*beta_s]);


Phi_0(1,1)=-2*1i*mu*beta_P*k_x;
Phi_0(1,2)=2*1i*mu*beta_P*k_x;
Phi_0(1,3)=1i*mu*(beta_s^2-k_x^2);
Phi_0(1,4)=1i*mu*(beta_s^2-k_x^2);


Phi_0(2,1)= beta_P;
Phi_0(2,2)=-beta_P;
Phi_0(2,3)= k_x;
Phi_0(2,4)= k_x;

Phi_0(3,1)=alpha_P;
Phi_0(3,2)=alpha_P;
Phi_0(3,3)=-alpha_s;
Phi_0(3,4)=alpha_s;


Phi_0(4,1)=k_x;
Phi_0(4,2)=k_x;
Phi_0(4,3)=-beta_s;
Phi_0(4,4)=beta_s;


% for iiii=1:4
%     [(-A)*Phi_0(:,iiii) V_0(iiii,iiii)*Phi_0(:,iiii) (-A)*Phi_0(:,iiii)-V_0(iiii,iiii)*Phi_0(:,iiii)]
% end
% pause


[a,indice]=sort(diag(real(V_0)),'descend');


for i_m=1:4
    Phi(:,i_m)=Phi_0(:,indice(4+1-i_m));
    lamda(i_m)=-V_0(indice(4+1-i_m),indice(4+1-i_m));
end

Psi=inv(Phi);
%MM=exp(lamda(1)*d)*Phi(:,1)*Psi(1,:)+exp(lamda(2)*d)*Phi(:,2)*Psi(2,:)+exp(lamda(3)*d)*Phi(:,3)*Psi(3,:)+exp(lamda(4)*d)*Phi(:,4)*Psi(4,:);


alpha_prime=Phi(:,2)*Psi(2,:)+Phi(:,3)*Psi(3,:)*exp((lamda(3)-lamda(2))*d)+Phi(:,4)*Psi(4,:)*exp((lamda(4)-lamda(2))*d);


Xi_prime=Psi(1:2,:)*Omega_moins;

%Omega_pm=MM*Omega_moins;

Omega_plus=[Phi(:,1) 0*Phi(:,2)]+alpha_prime*Omega_moins*inv(Xi_prime)*diag([exp((lamda(2)-lamda(1))*d) 1]);

T_back=inv(Xi_prime)*diag([exp(-lamda(1)*d) exp(-lamda(2)*d) ]);
% lamda
% lamda-lamda(1)
% pause



