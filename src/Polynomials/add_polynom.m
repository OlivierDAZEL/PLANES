function Q=add_polynom(P1,P2)


n=max([length(P1),length(P2)]);
Q=zeros(1,n);
Q(1:length(P1))=P1;
Q(1:length(P2))=Q(1:length(P2))+P2;
