model_data.lx=0.1;
model_data.ly=1;
model_data.nx=2;
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

% All the elements are DGM on TR6

elem.model=10*ones(nb.elements,1);
theta_DGM.nb=4;


