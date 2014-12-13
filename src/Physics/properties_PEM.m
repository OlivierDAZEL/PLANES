

% Biot densities with tortuosity effects
rho_12=-phi*air.rho*(alpha-1);
rho_11=rho_1-rho_12;
rho_2=phi*air.rho;
rho_22=rho_2-rho_12;

rho_22_til=phi^2*rho_eq_til;
rho_12_til=rho_2-rho_22_til;
rho_11_til=rho_1-rho_12_til;
rho_til=rho_11_til-((rho_12_til.^2)./rho_22_til);

gamma_til=phi*(rho_12_til./rho_22_til-(1-phi)/phi);
rho_s_til=rho_til+gamma_til.^2.*rho_eq_til;



% Biot in-vacuo elastic coefficients

N=(young*(1+1i*eta))/(2*(1+nu));
A_hat=(young*(1+1i*eta)*nu)/((1+nu)*(1-2*nu));
P_hat=A_hat+2*N;

% Biot 1956 elastic coefficients


R_til=K_eq_til*phi^2;
Q_til=((1-phi)/phi)*R_til;
P_til=P_hat+Q_til.^2./R_til;


% Parameters for energies

xi=(1-phi)/phi;
b_til=(1j*omega)*(rho_12-rho_12_til);
b_r=real(b_til);
rho_f_til=rho_til-rho_1;