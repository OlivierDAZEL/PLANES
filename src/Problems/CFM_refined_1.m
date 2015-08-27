
fid=fopen(name.file_input_FreeFem,'w');
fprintf(fid,'%s\n',name.file_msh);
fclose(fid);


% Call to FreeFem++ to create a msh File
system(['/usr/local/bin/FreeFem++ ' name.file_edp]);

%Importation of the msh File

[nb,nodes,elem,edge_msh]=msh_import(name.file_msh);

% All the elements are TR6
elem.model=1*ones(nb.elements,1);

label_boundary=25;
label_elem_ajoute=0;
model_elem_ajoute=11;
l_supp=6;
[nb,nodes,elem,edge_msh] = add_H_on_boundary_incompatible(l_supp,label_boundary,3,2,label_elem_ajoute,model_elem_ajoute,edge_msh,nodes,elem,nb);
theta_DGM.nb=16;
tilt=0*pi/theta_DGM.nb;