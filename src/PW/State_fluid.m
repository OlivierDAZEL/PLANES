function M=State_fluid(k_x,k_z,K)

k=sqrt(k_x^2+k_z^2);
M=[k_z -k_z;j*K*k^2 j*K*k^2];