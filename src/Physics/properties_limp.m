properties_jca

rho_2=phi*rho_0;

rho_22_til=phi^2*rho_eq_til;
rho_12_til=rho_2-rho_22_til;
rho_11_til=rho_1-rho_12_til;

gamma_til=phi*(rho_12_til./rho_22_til-(1-phi)/phi);

rho_til=rho_11_til-((rho_12_til.^2)./rho_22_til);
rho_s_til=rho_til+gamma_til.^2.*rho_eq_til;

rho_limp=rho_eq_til.*(rho_til./rho_s_til);
