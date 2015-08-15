model_data.lx=1.00e-2;
model_data.ly=5.00e-2;
model_data.nx=1;
model_data.ny=ceil(model_data.nx*model_data.ly/model_data.lx);



%solve.TR6=1;

fid=fopen(name.file_input_FreeFem,'w');
fprintf(fid,'%s\n',name.file_msh);
fprintf(fid,'%12.8f\n',model_data.lx);
fprintf(fid,'%12.8f\n',model_data.ly);
fprintf(fid,'%d\n',model_data.nx);
fprintf(fid,'%d\n',model_data.ny);
fclose(fid);

% Call to FreeFem++ to create a msh File
system(['/usr/local/bin/FreeFem++ ' name.file_edp]);


%Importation of the msh File

[nb,nodes,elem,edge_msh]=msh_import(name.file_msh);

% All the elements are TR6

elem.model=1*ones(nb.elements,1);

label_boundary=51
label_elem_ajoute=0
model_elem_ajoute=2
l_supp=0.05
[nb,nodes,elem,edge_msh] = add_H12_boundary(l_supp,label_boundary,label_elem_ajoute,model_elem_ajoute,edge_msh,nodes,elem,nb)
model_data.ly=model_data.ly+l_supp;


%aezeazezaeazaez