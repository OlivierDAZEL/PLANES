model_data.lx=1;
model_data.ly=1/2;

model_data.ny=ceil(model_data.nx*model_data.ly/model_data.lx);
theta_DGM.nb=16;


label_boundary=51;
model_data.label=[3 1 label_boundary 1];


[nb,nodes,elem,edge_msh]=createmshH12(model_data.lx,model_data.ly,model_data.nx,model_data.ny,model_data.label);
elem.model=2*ones(nb.elements,1);
% All the elements of the first domain are H12

label_elem_ajoute=0;
model_elem_ajoute=11;
l_supp=model_data.ly;
[nb,nodes,elem,edge_msh] = add_H12_boundary(l_supp,label_boundary,label_elem_ajoute,model_elem_ajoute,edge_msh,nodes,elem,nb);
model_data.ly=model_data.ly+l_supp;


