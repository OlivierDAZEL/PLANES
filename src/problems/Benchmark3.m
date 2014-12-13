theta=60*pi/180;

dyalu=1.0e-3;
nyalu=2;
dx=2e-3;
nx=ceil(dx*nyalu/dyalu);
dyporous=2.0e-2;
nyporous=ceil(dyporous*nyalu/dyalu);


labelexcitation=11;
labelalu=1001;
labelporous=5003;

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