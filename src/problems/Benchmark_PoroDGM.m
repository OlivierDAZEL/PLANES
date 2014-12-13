d=1;

fid=fopen(nom_fichier_input_FreeFem,'w');
fprintf(fid,'%12.8f\n',d);
fprintf(fid,'%d\n',n);
fprintf(fid,'%d\n',1000);
fprintf(fid,'%d\n',4001);
fclose(fid);

