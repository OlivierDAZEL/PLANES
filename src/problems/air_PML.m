d_x=0.1;
d_air=1;
d_PML=0.3;
nx=1;
n_air=ceil(nx*d_air/d_x);
n_PML=ceil(nx*d_PML/d_x);
labelPML=8010;


fid=fopen(name_file_input_FreeFem,'w');
fprintf(fid,'%s\n',name_file_msh);
fprintf(fid,'%12.8f\n',d_x);
fprintf(fid,'%12.8f\n',d_air);
fprintf(fid,'%12.8f\n',d_PML);
fprintf(fid,'%d\n',labelPML);
fprintf(fid,'%d\n',nx);
fprintf(fid,'%d\n',n_air);
fprintf(fid,'%d\n',n_PML);
fclose(fid);