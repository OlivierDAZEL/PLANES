theta=60*pi/180;

dyalu=1.0e-3;
nyalu=2;
dx=2e-2;
nx=10;
dyporous=2.0e-2;
nyporous=ceil(dyporous*nx/dx);


labelexcitation=13;
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


incident(1).typ=1001;
incident(1).thickness=1e-3;

transmission(1).typ=1001;
transmission(1).thickness=1e-3;

