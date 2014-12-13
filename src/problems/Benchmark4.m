theta=0*pi/180;

dyalu=1.0e-3;
nyalu=1;
dx=2e-3;
nx=ceil(dx*nyalu/dyalu);
dyporous=5.5e-2;
nyporous=ceil(dyporous*nyalu/dyalu);


labelexcitation=12;
labelalu=1001;
labelporous=5004;

fid=fopen(name_file_input_FreeFem,'w');
fprintf(fid,'%s\n',name_file_msh);
fprintf(fid,'%12.8f\n',dx);
fprintf(fid,'%d\n',nx);
fprintf(fid,'%12.8f\n',dyalu);
fprintf(fid,'%d\n',nyalu);
fprintf(fid,'%d\n',labelalu);
fprintf(fid,'%12.8f\n',dyporous);
fprintf(fid,'%d\n',nyporous);
fprintf(fid,'%d\n',labelporous);
fprintf(fid,'%d\n',labelexcitation);

fclose(fid);