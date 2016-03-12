
label_boundary=25;
data_model.label_boundaries=[3 1 label_boundary 1];
fid=fopen('FF.inp','w');
fprintf(fid,'%s\n',name.file_msh);
fprintf(fid,'%12.8f\n',data_model.lx);
fprintf(fid,'%12.8f\n',data_model.ly/2);
fprintf(fid,'%d\n',data_model.nx);
fprintf(fid,'%d\n',data_model.ny);
fprintf(fid,'%d\n',data_model.label_boundaries(1));
fprintf(fid,'%d\n',data_model.label_boundaries(2));
fprintf(fid,'%d\n',data_model.label_boundaries(3));
fprintf(fid,'%d\n',data_model.label_boundaries(4));
fclose(fid);

% Call to FreeFem++ to create a msh File
system(['/usr/local/bin/FreeFem++ ' name.file_edp]);

%Importation of the msh File

[nb,nodes,elem,edge_msh]=msh_import(name.file_msh);

% All the elements are TR6
elem.model=1*ones(nb.elements,1);


label_elem_ajoute=0;
model_elem_ajoute=11;
l_supp=data_model.ly/2;
[nb,nodes,elem,edge_msh] = add_H12_boundary(l_supp,label_boundary,label_elem_ajoute,model_elem_ajoute,edge_msh,nodes,elem,nb);
%data_model.ly=data_model.ly+l_supp;