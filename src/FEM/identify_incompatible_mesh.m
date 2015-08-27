

for ie=1:size(edges.incompatible,1)
    if (nodes(edges.incompatible(ie,1),1))>(nodes(edges.incompatible(ie,2),1))
        temp=edges.incompatible(ie,1);
        edges.incompatible(ie,1)=edges.incompatible(ie,2);
        edges.incompatible(ie,2)=temp;
    end
end


% For each edge in internal_FEMFEM x(n1)<x(n2)


% Finds the coordinates of all the nodes of the incompatible

nodes_incompatible=[edges.incompatible(:,1);edges.incompatible(:,2)];






[~,index_condensed,~]=unique(nodes(nodes_incompatible,:),'rows');

nodes_incompatible=nodes_incompatible(index_condensed);


% nodes incompatible are ordered by increasing x
x_nodes_incompatible=nodes(nodes_incompatible,1);


[~,i_kept_incompatible,i_incompatible]=unique(x_nodes_incompatible);


nodes_incompatible_edge=nodes_incompatible(i_kept_incompatible);


corresp=[];

for ii=1:length(nodes_incompatible)
    corresp(nodes_incompatible(ii))=nodes_incompatible_edge(i_incompatible(ii));
end




nb_nodes_condensed=length(nodes_incompatible_edge);
edges.flux=zeros(nb_nodes_condensed-1,4);
edges.flux(1:nb_nodes_condensed-1,1)=nodes_incompatible_edge(1:end-1)';
edges.flux(1:nb_nodes_condensed-1,2)=nodes_incompatible_edge(2:end)';



x_nodes_edge=nodes(nodes_incompatible_edge,1);




for ie=1:size(edges.incompatible,1)
   x_min=nodes(edges.incompatible(ie,1),1);
   x_max=nodes(edges.incompatible(ie,2),1);
   i_nodes_condensed_min=find(x_nodes_edge==x_min);
   i_nodes_condensed_max=find(x_nodes_edge==x_max);
   for i_temp=i_nodes_condensed_min:i_nodes_condensed_max-1
      if  edges.flux(i_temp,3)==0
          edges.flux(i_temp,3)=edges.incompatible(ie,3);
      else
          edges.flux(i_temp,4)=edges.incompatible(ie,3);
   end
   
   end
end

nb.flux=size(edges.flux,1);


edges=rmfield(edges,'incompatible');
