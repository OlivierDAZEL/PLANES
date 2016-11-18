%solve.TR6=1;
label_boundary=25;
data_model.label=[3 1 label_boundary 1];


[nb,nodes,elem,edge_msh]=createmshH12(data_model.lx,data_model.ly/2,data_model.nx,data_model.ny,data_model.label);

% All the elements are H12

elem.model=11*ones(nb.elements,1);

label_elem_ajoute=8010;
model_elem_ajoute=11;
l_supp=data_model.ly/2;
[nb,nodes,elem,edge_msh] = add_H12_boundary(l_supp,label_boundary,label_elem_ajoute,model_elem_ajoute,edge_msh,nodes,elem,nb);



