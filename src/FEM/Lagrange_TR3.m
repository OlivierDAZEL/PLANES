function Poly = Lagrange_TR3(nodes,num)

A=[ones(3,1) nodes];
b=zeros(3,1);
b(num)=1;
X=A\b;
Poly=[X(1) X(2);X(3) 0];

end

