fid=fopen('FF.inp','w');
fprintf(fid,'%s\n',name.file_msh);
fprintf(fid,'%12.8f\n',data_model.lx);
fprintf(fid,'%12.8f\n',data_model.ly);
fprintf(fid,'%d\n',data_model.nx);
fprintf(fid,'%d\n',data_model.ny);
fclose(fid);

% Call to FreeFem++ to create a msh File
system(['/usr/local/bin/FreeFem++ ' name.file_edp]);

%Importation of the msh File

[nb,nodes,elem,edge_msh]=msh_import(name.file_msh);

% All the elements are TR6
elem.model=10*ones(nb.elements,1);

