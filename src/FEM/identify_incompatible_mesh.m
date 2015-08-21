

for ie=1:size(edges.internal_FEMFEM,1)
    if (nodes(edges.internal_FEMFEM(ie,1),1))>(nodes(edges.internal_FEMFEM(ie,2),1))
        temp=edges.internal_FEMFEM(ie,1);
        edges.internal_FEMFEM(ie,1)=edges.internal_FEMFEM(ie,2);
        edges.internal_FEMFEM(ie,2)=temp;
    end
end
% For each edge in internal_FEMFEM x(n1)<x(n2)


% Finds the coordinates of all the nodes of the incompatible

nodes_incompatible=[edges.internal_FEMFEM(:,1);edges.internal_FEMFEM(:,2)];



[~,index_condensed,~]=unique(nodes(nodes_incompatible,:),'rows');

nodes_incompatible_condensed=nodes_incompatible(index_condensed);
[~,i_temp]=sort(nodes(nodes_incompatible_condensed,1));
nodes_incompatible_condensed=nodes_incompatible_condensed(i_temp);
% nodes incompatible are ordered by increasing x
x_nodes_incompatible_condensed=nodes(nodes_incompatible_condensed,1);
nb_nodes_condensed=length(nodes_incompatible_condensed);


edges.flux=zeros(nb_nodes_condensed-1,4);
edges.flux(1:nb_nodes_condensed-1,1)=nodes_incompatible_condensed(1:end-1);
edges.flux(1:nb_nodes_condensed-1,2)=nodes_incompatible_condensed(2:end);

for ie=1:size(edges.internal_FEMFEM,1)
   x_min=nodes(edges.internal_FEMFEM(ie,1),1);
   x_max=nodes(edges.internal_FEMFEM(ie,2),1);
   i_nodes_condensed_min=find(x_nodes_incompatible_condensed==x_min);
   i_nodes_condensed_max=find(x_nodes_incompatible_condensed==x_max);
   for i_temp=i_nodes_condensed_min:i_nodes_condensed_max-1
      if  edges.flux(i_temp,3)==0
          edges.flux(i_temp,3)=edges.internal_FEMFEM(ie,3);
      else
          edges.flux(i_temp,4)=edges.internal_FEMFEM(ie,3);
   end
   
   end
end

nb.flux=size(edges.flux,1);

