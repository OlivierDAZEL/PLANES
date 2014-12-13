
x=[];
y=[];
% figure(567)
% hold on
% for ie=1:nb_elements
% 
%     c=X_EF_fortran(3*kconec(ie,:));
%     vertices=[vcor(kconec(ie,:),:)'];    
%     plot(vertices(2,:),abs(c),'.');
%      
% end
% title('Pression EF (module)')


% hold on
% xx=linspace(0,1,200);
% plot(xx,2*abs(cos(omega*(xx-1)/c_0)),'r');





figure(5679)
hold on
for ie=1:nb_elements

    c=X_EF_fortran(3*elements(ie,1:6));
    vertices=[nodes(elements(ie,1:6),:)'];    
    plot(vertices(2,:),abs(c),'.');
     
end
title('Pression EF (module)')




figure(5678)
hold on
for ie=1:nb_elements

    c=X_EF_fortran(3*(elements(ie,1:6)-1)+2);
    vertices=[nodes(elements(ie,1:6),:)'];    
    plot(vertices(2,:),abs(c),'.');
     
end
title('u_y EF (module)')












% figure(56780)
% hold on
% for ie=1:nb_elements
% 
%     c=X_EF_fortran(3*(kconec(ie,:)-1)+2);
%     vertices=[vcor(kconec(ie,:),:)'];    
%     plot(vertices(2,:),angle(c),'.');
%      
% end
% title('u_y EF (phase)')
% 
% 
% figure(56790)
% hold on
% for ie=1:nb_elements
% 
%     c=X_EF_fortran(3*kconec(ie,:));
%     vertices=[vcor(kconec(ie,:),:)'];    
%     plot(vertices(2,:),angle(c),'.');
%      
% end
% title('Pression EF (phase)')
