theta=60*pi/180;
dxair=2e-2;
dyair=5.5e-2;
nx=10;
ny=10;
labelexcitation=12;
labelporous=5004;
r=2.5e-3;
labelelastic=1001;
ninclusion=2*nx;


fid=fopen(name_file_input_FreeFem,'w');
fprintf(fid,'%s\n',name_file_msh);
fprintf(fid,'%12.8f\n',dxair);
fprintf(fid,'%12.8f\n',dyair);
fprintf(fid,'%d\n',nx);
fprintf(fid,'%d\n',ny);
fprintf(fid,'%d\n',labelexcitation);
fprintf(fid,'%d\n',labelporous);
fprintf(fid,'%12.8f\n',r);
fprintf(fid,'%d\n',labelelastic);
fprintf(fid,'%d\n',ninclusion);
fclose(fid);