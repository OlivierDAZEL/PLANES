%disp('Applying periodicity condtions')



edge_left= find(periodicity(:,4)==98);
edge_right=find(periodicity(:,4)==99);

node_left=unique([periodicity(edge_left,1);periodicity(edge_left,2);periodicity(edge_left,6)]);
[temp,i_left]=sort(nodes(node_left,2));
node_left=node_left(i_left);

node_right=unique([periodicity(edge_right,1);periodicity(edge_right,2);periodicity(edge_right,6)]);
[temp,i_right]=sort(nodes(node_right,2));
node_right=node_right(i_right);

delta=exp(-1i*k_x*period);

dof_left= dof_A([3*(node_left -1)+1;3*(node_left -1)+2;3*(node_left -1)+3]);
dof_right=dof_A([3*(node_right-1)+1;3*(node_right-1)+2;3*(node_right-1)+3]);

dof_left= dof_left (find(dof_left));
dof_right=dof_right(find(dof_right));


for ii=1:length(dof_left)
    
        A(:,dof_left(ii))=A(:,dof_left(ii))+delta*A(:,dof_right(ii));
        A(:,dof_right(ii))=0;

        
        A(dof_left(ii),:)=A(dof_left(ii),:)+A(dof_right(ii),:)/delta;
        F(dof_left(ii))=F(dof_left(ii))+F(dof_right(ii))/delta;

        A(dof_right(ii),:)=0;
        
        F(dof_right(ii))=0;
        A(dof_right(ii),dof_left(ii))=delta;
        A(dof_right(ii),dof_right(ii))=-1;
        
    
end