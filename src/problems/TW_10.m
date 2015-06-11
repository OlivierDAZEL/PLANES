period=20e-3;
thicknessporous=period;
labelporous=5005;
labelinclusion=1003;
radius=5.0e-3;
name_file_TW_export='FEM_Abs_d20_angpi90_type4alr5.dat'


nx=12;
ny=nx;
ninclusion=ceil(nx*2*pi*radius/period);


fid=fopen(name_file_input_FreeFem,'w');
fprintf(fid,'%s\n',name_file_msh);
fprintf(fid,'%12.8f\n',period);
fprintf(fid,'%12.8f\n',thicknessporous);
fprintf(fid,'%d\n',labelporous);
fprintf(fid,'%d\n',labelinclusion);
fprintf(fid,'%12.8f\n',radius);
fprintf(fid,'%d\n',nx);
fprintf(fid,'%d\n',ny);
fprintf(fid,'%d\n',ninclusion);
fclose(fid);


