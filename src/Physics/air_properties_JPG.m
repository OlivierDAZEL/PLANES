



%% Constants
air.T=20;                                     %[?C]
air.P=1.01325e5;                              %[Pa]
air.gamma=1.400;                                %[1]
air.mmol=.29e-1;                                %[kg.mol^-1]
air.mu=.1839e-4; % Dynamic viscosity           %(kg/m/s)
air.NPR=0.710;
air.lambda=0.0262;                              %[W.m^-1.K^-1]

%% Air properties
air.rho=1.213;                                %[kg.m^-3]
air.c=341.9730829285405;                              %[m.s^-1]

air.K=air.c^2*air.rho;
air.Z=air.rho*air.c;

air.C_p=1006;
air.C_v=air.C_p/air.gamma;

air.nu=air.mu/air.rho;
air.nu_prime=air.nu/air.NPR;
