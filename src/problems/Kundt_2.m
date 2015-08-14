model_data.lx=1.00e-2;
model_data.ly=5.00e-2;
model_data.nx=1;
model_data.ny=ceil(model_data.nx*model_data.ly/model_data.lx);
mdel_data.label=[60 1 1 1];


[nb,nodes,elem,edge_msh]=createmshH12(model_data.lx,model_data.ly,model_data.nx,model_data.ny,mdel_data.label);

% All the elements are H12

elem.model=2*ones(nb.elements,1);




%aezeazezaeazaez