d_x=0.1;
d_air=0.8;
d_EF=0.8;
nx=1;
n_air=ceil(nx*d_air/d_x);
n_EF=ceil(nx*d_EF/d_x);
labelEF=2004;

fid=fopen(name_file_input_FreeFem,'w');
fprintf(fid,'%s\n',name_file_msh);
fprintf(fid,'%12.8f\n',d_x);
fprintf(fid,'%12.8f\n',d_air);
fprintf(fid,'%12.8f\n',d_EF);
fprintf(fid,'%d\n',labelEF);
fprintf(fid,'%d\n',nx);
fprintf(fid,'%d\n',n_air);
fprintf(fid,'%d\n',n_EF);
fclose(fid);