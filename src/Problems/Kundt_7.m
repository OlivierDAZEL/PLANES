model_data.lx=0.1;
model_data.ly=1;
model_data.nx=2;
model_data.ny=ceil(model_data.nx*model_data.ly/model_data.lx);
label_boundary=25;
model_data.label_boundaries=[3 1 label_boundary 1];
theta_DGM.nb=4;

fid=fopen(name.file_input_FreeFem,'w');
fprintf(fid,'%s\n',name.file_msh);
fprintf(fid,'%12.8f\n',model_data.lx);
fprintf(fid,'%12.8f\n',model_data.ly);
fprintf(fid,'%d\n',model_data.nx);
fprintf(fid,'%d\n',model_data.ny);
fprintf(fid,'%d\n',model_data.label_boundaries(1));
fprintf(fid,'%d\n',model_data.label_boundaries(2));
fprintf(fid,'%d\n',model_data.label_boundaries(3));
fprintf(fid,'%d\n',model_data.label_boundaries(4));
fclose(fid);

% Call to FreeFem++ to create a msh File
system(['/usr/local/bin/FreeFem++ ' name.file_edp]);


%Importation of the msh File

[nb,nodes,elem,edge_msh]=msh_import(name.file_msh);


% All the elements are TR6

elem.model=1*ones(nb.elements,1);


label_elem_ajoute=0;
model_elem_ajoute=11;
l_supp=model_data.ly;
[nb,nodes,elem,edge_msh] = add_H12_boundary(l_supp,label_boundary,label_elem_ajoute,model_elem_ajoute,edge_msh,nodes,elem,nb);
model_data.ly=model_data.ly+l_supp;


%aezeazezaeazaez