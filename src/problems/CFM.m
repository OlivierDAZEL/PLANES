lx=1.00e-2;
lyFEM=10.00e-2;

lyDGM=lyFEM;
nx=1;
nyFEM=ceil(nx*lyFEM/lx);
nyDGM=ceil(nx*lyDGM/lx);

% fid=fopen(name_file_input_FreeFem,'w');
% fprintf(fid,'%s\n',name_file_msh);
% fprintf(fid,'%12.8f\n',lx);
% fprintf(fid,'%12.8f\n',ly);
% fprintf(fid,'%d\n',nx);
% fprintf(fid,'%d\n',ny);
% fclose(fid);


