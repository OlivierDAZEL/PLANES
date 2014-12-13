
dyalu=1.0e-3;

nyalu=8;

dx=1e-3;
nx=ceil(dx*nyalu/dyalu);
%dyalu*(nx/nyalu);




labelexcitation=20045;
labelelement=1001;

fid=fopen(nom_fichier_input_FreeFem,'w');
fprintf(fid,'%s\n',nom_fichier_msh);
fprintf(fid,'%12.8f\n',dx);
fprintf(fid,'%12.8f\n',dyalu);
fprintf(fid,'%d\n',nx);
fprintf(fid,'%d\n',nyalu);
fprintf(fid,'%d\n',labelexcitation);
fprintf(fid,'%d\n',labelelement);
fclose(fid);