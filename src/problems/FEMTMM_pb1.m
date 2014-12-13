theta=80*pi/180;


delta_x_MMT=0;
delta_y_MMT=1e-2;
dx=10e-2;
dybas=10e-2;
dyhaut=1e-2;

nx=20;
nybas=20;
nyhaut=2;

nyporous=nyhaut*ceil(delta_y_MMT/dyhaut);

labelexcitation=10;
labeltermination=20;
termination=1;
label_porous=4003;

% labelexcitation=10;
% labelelement=4003;

fid=fopen(name_file_input_FreeFem,'w');
fprintf(fid,'%s\n',name_file_msh);
fprintf(fid,'%12.8f\n',dx);
fprintf(fid,'%12.8f\n',dybas);
fprintf(fid,'%12.8f\n',dyhaut);
fprintf(fid,'%d\n',nx);
fprintf(fid,'%d\n',nybas);
fprintf(fid,'%d\n',nyhaut);
fprintf(fid,'%d\n',nyporous);
fprintf(fid,'%d\n',labelexcitation);
fprintf(fid,'%d\n',labeltermination);
fprintf(fid,'%d\n',label_porous);
fprintf(fid,'%12.8f\n',delta_y_MMT);
fclose(fid);

% Including the air incident medium (layer #1)
nb_layers=3;
multilayer(1).d=dybas;
multilayer(1).mat=0;
multilayer(2).d=delta_y_MMT;
multilayer(2).mat=label_porous;
multilayer(3).d=dyhaut;
multilayer(3).mat=0;
% Termination condition // 0 for rigid backing 1 for radiation
%termination=0;

% Number of waves (two ways) in each layer
% For the resolution, the incident waves is included in the system and put to RHS at the
% end of the procedure

% Addition of a new layer for the incident medium
l0.d=0;
l0.mat=0;
multilayer=[l0 multilayer];
nb_layers=nb_layers+1;

compute_number_PW_TMM
