lx=1.00e-2;
ly=1.00e-1;
nx=3;
ny=ceil(nx*ly/lx);


fid=fopen(name_file_input_FreeFem,'w');
fprintf(fid,'%s\n',name_file_msh);
fprintf(fid,'%12.8f\n',lx);
fprintf(fid,'%12.8f\n',ly);
fprintf(fid,'%d\n',nx);
fprintf(fid,'%d\n',ny);
fclose(fid);


