function [vh,vq] = TR6_fluid_axi(coord_e)



% Gauss Points Data
%------------------------------

npg=6;
a=0.445948490915965;
b=0.091576213509771;
P1=0.111690794839005;
P2=0.054975871827661;


ksig=[  a, a;
    1-2*a, a;
    a,   1-2*a;
    b, b;
    1-2*b, b;
    b, 1-2*b];
w_g=[P1,P1,P1,P2,P2,P2];


% Matrices initialization
%-----------------------

vh=zeros(6,6);
vq=zeros(6,6);



% boucle sur les PG
%------------------
for ipg=1:npg,
    
    ksi=ksig(ipg,1);
    eta=ksig(ipg,2);
    lambda=1-ksi-eta;
    
    Phi =  [ -lambda*(1-2*lambda)  4*ksi*lambda  -ksi*(1-2*ksi) 4*ksi*eta -eta*(1-2*eta) 4*eta*lambda];
    dPhi  = [ 1-4*lambda 4*(lambda-ksi) -1+4*ksi 4*eta 0 -4*eta; ...
        1-4*lambda -4*ksi 0 4*ksi -1+4*eta 4*(lambda-eta)];
    J = dPhi*coord_e';
    
    r=Phi*coord_e(1,:)';
    
    weight = r*w_g(ipg) * det(J);
    
    vB =J\dPhi;
    
    Gd=[vB(1,1) vB(1,2) vB(1,3) vB(1,4) vB(1,5) vB(1,6);
        vB(2,1) vB(2,2) vB(2,3) vB(2,4) vB(2,5) vB(2,6)];
    
    vh=vh+Gd'*Gd*weight;
    vq=vq+Phi'*Phi*weight;
    
end;   % fin de boucle sur les points de gauss.



