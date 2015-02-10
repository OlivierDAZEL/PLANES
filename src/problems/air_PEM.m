d_x=0.1;
d_air=1;
d_PEM=1;
nx=3;
n_air=ceil(nx*d_air/d_x);
n_PEM=ceil(nx*d_PEM/d_x);
labelPEM=5004;


fid=fopen(name_file_input_FreeFem,'w');
fprintf(fid,'%s\n',name_file_msh);
fprintf(fid,'%12.8f\n',d_x);
fprintf(fid,'%12.8f\n',d_air);
fprintf(fid,'%12.8f\n',d_PEM);
fprintf(fid,'%d\n',labelPEM);
fprintf(fid,'%d\n',nx);
fprintf(fid,'%d\n',n_air);
fprintf(fid,'%d\n',n_PEM);
fclose(fid);