lx=0.2;
ly=0.2;
ax=0.2;
ay=0.2;
bx=0.2;
by=0.2;
rayon=0.05;
nlx=2;
nly=2;
nax=2;
nay=2;
nbx=2;
nby=2;
nsrc=16;



fid=fopen(nom_fichier_input_FreeFem,'w');
fprintf(fid,'%s\n',nom_fichier_msh);
fprintf(fid,'%12.8f\n',lx);
fprintf(fid,'%12.8f\n',ly);
fprintf(fid,'%12.8f\n',ax);
fprintf(fid,'%12.8f\n',ay);
fprintf(fid,'%12.8f\n',bx);
fprintf(fid,'%12.8f\n',by);
fprintf(fid,'%12.8f\n',rayon);
fprintf(fid,'%d\n',nlx);
fprintf(fid,'%d\n',nly);
fprintf(fid,'%d\n',nax);
fprintf(fid,'%d\n',nay);
fprintf(fid,'%d\n',nbx);
fprintf(fid,'%d\n',nby);
fprintf(fid,'%d\n',nsrc);
fclose(fid);

