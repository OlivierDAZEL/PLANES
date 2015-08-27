
fid=fopen(name.file_input_FreeFem,'w');
fprintf(fid,'%s\n',name.file_msh);
fclose(fid);


% Call to FreeFem++ to create a msh File
system(['/usr/local/bin/FreeFem++ ' name.file_edp]);

%Importation of the msh File

[nb,nodes,elem,edge_msh]=msh_import(name.file_msh);

% All the elements are TR6
elem.model=1*ones(nb.elements,1);
