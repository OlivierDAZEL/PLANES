dyalu=1.2e-3;
nyalu=10;
dyair=1e-2;
dx=1.2e-3;

nx=ceil(nyalu*dx/dyalu)
nyair=ceil(nyalu*dyair/dyalu);


labelexcitation=10000;
labelelement=1001;


fid=fopen(nom_fichier_input_FreeFem,'w');
fprintf(fid,'%s\n',nom_fichier_msh);
fprintf(fid,'%12.8f\n',dx);
fprintf(fid,'%12.8f\n',dyalu);
fprintf(fid,'%12.8f\n',dyair);
fprintf(fid,'%d\n',nx);
fprintf(fid,'%d\n',nyalu);
fprintf(fid,'%d\n',nyair);
fprintf(fid,'%d\n',labelexcitation);
fprintf(fid,'%d\n',labelelement);
fclose(fid);