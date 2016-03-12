data_model.lx=1;
data_model.ly=1;
data_model.nx=1;
data_model.ny=ceil(data_model.nx*data_model.ly/data_model.lx);
data_model.theta_DGM.nb=4;
data_model.tilt=0*pi/data_model.theta_DGM.nb;


%solve.TR6=1;

fid=fopen(name.file_input_FreeFem,'w');
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

% All the elements are DGM on TR

elem.model=10*ones(nb.elements,1);


