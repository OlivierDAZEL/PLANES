function f=int_edge_2vectorielle(k_1,k_2,a,b,centres)

% return a matrix
% f(i,j)=\int_a^b e^{(k_i*(x-centres(1)+k_j*(x-centres(2)))}

a_1=a-centres(1:2,1);
a_2=a-centres(1:2,2);

k_1
k_2

t=(b-a)/norm(b-a)

k1t=k_1.'*t
k2t=k_2.'*t

temp=(k1t*ones(1,length(k_2))+ones(length(k_1),1)*k2t.')

[i_temp,j_temp]=find(abs(temp)<1e-6);
temp(i_temp,j_temp)=0;

temp0=temp.*temp==0
tempn0=temp~=0

temp=temp+temp0; % Pour la division

norm(b-a)
((exp(norm(b-a)*temp)-1)./temp)

f=norm(b-a)*temp0+tempn0.*((exp(norm(b-a)*temp)-1)./temp)

%%%%
%size(f);
%size(exp((k_1.'*a_1)*ones(1,length(k_2))));
%ones(length(k_1),1);
%(k_2.'*a_2);


f=f.*exp((k_1.'*a_1)*ones(1,length(k_2))).*exp(ones(length(k_1),1)*(k_2.'*a_2).');

