function M=Phi_Biot(nx,ny,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega)


M=[nx;ny ;mu_1*nx;mu_1*ny ;-(delta_1/omega)*(A_hat+2*N*nx.^2);-(delta_1/omega)*(2*N*nx.*ny)    ;-(delta_1/omega)*(A_hat+2*N*ny.^2);-(delta_1/omega)*(K_eq_til*mu_1)*ones(1,length(nx))];
M=[M [nx;ny ;mu_2*nx;mu_2*ny ;-(delta_2/omega)*(A_hat+2*N*nx.^2);-(delta_2/omega)*(2*N*nx.*ny)    ;-(delta_2/omega)*(A_hat+2*N*ny.^2);-(delta_2/omega)*(K_eq_til*mu_2)*ones(1,length(nx))]];
M=[M [ny;-nx;mu_3*ny;-mu_3*nx;-(delta_3/omega)*(2*N*nx.*ny)     ; (delta_3/omega)*(N*(nx.^2-ny.^2)); (delta_3/omega)*(2*N*nx.*ny)                         ;0*ones(1,length(nx))]];




% M(3,4)=ny;
% M(4,4)=-nx;
% 
% M(5,5)=ny^2;
% M(6,5)=-nx*ny;
% M(7,5)=nx^2;

M5=(M(5,:)+M(7,:))/2;
M7=(M(5,:)-M(7,:))/2;

M(5,:)=M5;
M(7,:)=M7;