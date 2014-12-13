period=50e-2;
thicknessplate=1e-3;
thicknessporous=2e-2;
labelplate=1001;
labelporous=5003;
hbottom=5e-2;
htop=5e-2;
ampsinus=1*thicknessporous;
ampx2=1;
nboscill=1;

nplate=1;


fid=fopen(name_file_input_FreeFem,'w');
fprintf(fid,'%s\n',name_file_msh);
fprintf(fid,'%12.8f\n',period);
fprintf(fid,'%12.8f\n',thicknessplate);
fprintf(fid,'%12.8f\n',thicknessporous);
fprintf(fid,'%d\n',labelplate);
fprintf(fid,'%d\n',labelporous);
fprintf(fid,'%12.8f\n',hbottom);
fprintf(fid,'%12.8f\n',htop);
fprintf(fid,'%12.8f\n',ampsinus);
fprintf(fid,'%12.8f\n',ampx2);
fprintf(fid,'%12.8f\n',nboscill);
fprintf(fid,'%d\n',nplate);
fclose(fid);


delta_y_MMT=thicknessporous;
delta_x_MMT=0;

multilayer_femtmm(1).d=thicknessporous;
multilayer_femtmm(1).mat=labelporous;



nb_layers=3;
multilayer(1).d=thicknessplate;
multilayer(1).mat=labelplate;
multilayer(2).d=thicknessporous;
multilayer(2).mat=labelporous;
multilayer(3).d=thicknessplate;
multilayer(3).mat=labelplate;
% Termination condition // 0 for rigid backing 1 for radiation
termination=1;

% Number of waves (two ways) in each layer
% For the resolution, the incident waves is included in the system and put to RHS at the
% end of the procedure

% Addition of a new layer for the incident medium
l0.d=0;
l0.mat=0;
multilayer=[l0 multilayer];
nb_layers=nb_layers+1;

compute_number_PW_TMM