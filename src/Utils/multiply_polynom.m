function Q=multiply_polynom(P1,P2)

Q=zeros(1,(length(P1)-1)+(length(P2)-1)+1);
for ii=1:length(P1)
   for jj=1:length(P2)
      Q(ii+jj-1)=Q(ii+jj-1)+P1(ii)*P2(jj); 
   end
end

