k_air=omega/air.c;
k_x=k_air*sin(theta);
k_z=k_air*cos(theta);

nb_Bloch_waves=ceil((period/(2*pi))*(3*real(k_air)-k_x))+10;

nb_Bloch_waves=floor((period/(2*pi))*(3*real(k_air)-k_x))+5;

%nb_Bloch_waves=0
if nb_R~=0
    nb_R=2*nb_Bloch_waves+1;
end
if nb_T~=0
    nb_T=2*nb_Bloch_waves+1;
end

temp=[];
temp(1:2:2*nb_Bloch_waves)=1:nb_Bloch_waves;
temp(2:2:2*nb_Bloch_waves+1)=-(1:nb_Bloch_waves);
temp=[0 temp];

vec_k_x=k_x+temp*(2*pi/period);
vec_k_x_t=k_x+temp*(2*pi/period);

vec_k_z=sqrt(k_air^2-vec_k_x.^2);
vec_k_z=real(vec_k_z)-1i*imag(vec_k_z);
vec_k_z_t=sqrt(k_air^2-vec_k_x_t.^2);
vec_k_z_t=real(vec_k_z_t)-1i*imag(vec_k_z_t);
