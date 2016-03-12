%solve.TR6=1;

data_model.label=[3 1 1 1];


[nb,nodes,elem,edge_msh]=createmshH12(data_model.lx,data_model.ly,data_model.nx,data_model.ny,data_model.label);

% All the elements are H12

elem.model=2*ones(nb.elements,1);



