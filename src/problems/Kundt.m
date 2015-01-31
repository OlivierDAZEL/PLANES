lx=1.00;
ly=1e-2;
nx=30;
ny=ceil(nx*ly/lx);


fid=fopen(name_file_input_FreeFem,'w');
fprintf(fid,'%s\n',name_file_msh);
fprintf(fid,'%12.8f\n',lx);
fprintf(fid,'%12.8f\n',ly);
fprintf(fid,'%d\n',nx);
fprintf(fid,'%d\n',ny);
fclose(fid);


