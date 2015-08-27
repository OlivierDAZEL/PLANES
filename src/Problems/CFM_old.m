lAB=1.00e0;
lBC=1.00e0;
lCDhori=1.00e0;
lCDvert=1.00e0;
lDE=1.00e0;
lFE=lAB+lCDhori;
lGFvert=1.00e0;
lHG=1.00e0;
lAH=lBC+lCDvert+lDE-lGFvert;

nAB=5;
nBC=nAB*ceil(lBC/lAB);
nCD=nAB*ceil(sqrt(lCDhori^2+lCDvert^2)/lAB);
nDE=nAB*ceil(lDE/lAB);
nFE=nAB*ceil(lFE/lAB);
nGF=nAB*ceil(sqrt(lGFvert^2+lAB^2)/lAB);
nHG=nAB*ceil(lHG/lAB);
nAH=nAB*ceil(lHG/lAB);



fid=fopen(name.file_input_FreeFem,'w');
fprintf(fid,'%s\n',name.file_msh);
fprintf(fid,'%12.8f\n',lAB);
fprintf(fid,'%12.8f\n',lBC);
fprintf(fid,'%12.8f\n',lCDhori);
fprintf(fid,'%12.8f\n',lCDvert);
fprintf(fid,'%12.8f\n',lDE);
fprintf(fid,'%12.8f\n',lFE);
fprintf(fid,'%12.8f\n',lGFvert);
fprintf(fid,'%12.8f\n',lHG);

fprintf(fid,'%d\n',nAB);
fprintf(fid,'%d\n',nBC);
fprintf(fid,'%d\n',nCD);
fprintf(fid,'%d\n',nDE);
fprintf(fid,'%d\n',nFE);
fprintf(fid,'%d\n',nGF);
fprintf(fid,'%d\n',nHG);
fprintf(fid,'%d\n',nAH);
fclose(fid);

% Call to FreeFem++ to create a msh File
system(['/usr/local/bin/FreeFem++ ' name.file_edp]);


%Importation of the msh File

[nb,nodes,elem,edge_msh]=msh_import(name.file_msh);

% All the elements are TR6
elem.model=1*ones(nb.elements,1);
