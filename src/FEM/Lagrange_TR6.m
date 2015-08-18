function Poly = Lagrange_TR6(nodes,num)

A=[ones(6,1) nodes nodes(:,1).^2 nodes(:,2).^2 nodes(:,1).*nodes(:,2)];
% X(1)+X(2)*x+X(3)*y+X(4)*x^2+X(5)*y^2+X(6)*xy
b=zeros(6,1);
b(num)=1;
X=A\b;
Poly=[X(1) X(2) X(4);X(3) X(6) 0;X(5) 0 0];

end

