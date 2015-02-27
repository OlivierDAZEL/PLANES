function f=int_edge_1vectorielle(k_1,a,b,centres)

% return a matrix
% f(i)=\int_a^b e^{(k_i*(x-centres(1)))}

a_1=a-centres(1:2,1);

t=(b-a)/norm(b-a);

k1t=k_1.'*t;

temp=(k1t.');

[i_temp,j_temp]=find(abs(temp)<1e-6);
temp(i_temp,j_temp)=0;

temp0=temp.*temp==0;
tempn0=temp~=0;

temp=temp+temp0; % Pour la division

f=norm(b-a)*temp0+tempn0.*((exp(norm(b-a)*temp)-1)./temp);

%%%%
%size(f);
%size(exp((k_1.'*a_1)*ones(1,length(k_2))));
%ones(length(k_1),1);
%(k_2.'*a_2);



f=(f.').*exp((k_1.'*a_1));
