

%%  Johnson et al model for rho_eq_til

omega_0=sig*phi/(air.rho*alpha);
omega_infty=(sig*phi*LCV)^2/(4*air.mu*air.rho*alpha^2);
F_JKD=sqrt(1+1j*omega/omega_infty);
rho_eq_til=(air.rho*alpha/phi)*(1+(omega_0/(1j*omega))*F_JKD);
alpha_til=phi*rho_eq_til./air.rho;


%%  Champoux-Allard model for K_eq_til

omega_prime_infty=(16*air.nu_prime)/(LCT^2);
F_prime_CA=sqrt(1+1j*omega./omega_prime_infty);
alpha_prime_til=1+omega_prime_infty.*F_prime_CA./(2*1j*omega);
K_eq_til=(air.gamma*air.P./phi)./(air.gamma-(air.gamma-1)./alpha_prime_til);



