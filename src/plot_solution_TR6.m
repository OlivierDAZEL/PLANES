

for ie=1:nb_elements
     coord_e=nodes(elements(ie,:),:);
     sol_e=sol(elements(ie,:));
     xdata=[coord_e( [1 2 6],1) coord_e( [2 3 4],1) coord_e( [2 4 6],1) coord_e( [4 5 6],1)];
     ydata=[coord_e( [1 2 6],2) coord_e( [2 3 4],2) coord_e( [2 4 6],2) coord_e( [4 5 6],2)];
     zdata=abs([sol_e( [1 2 6]) sol_e( [2 3 4]) sol_e( [2 4 6]) sol_e( [4 5 6])]);
     patch(xdata,ydata,zdata,'EdgeColor','none')
     line([coord_e([1 3 5 1],1)],[coord_e([1 3 5 1],2)],'Color','k')
end
     

