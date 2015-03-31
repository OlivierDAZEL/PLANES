lx=10.00e-2;
ly=10.00e-2;
nx=6;

ny=ceil(nx*ly/lx);


fid=fopen(name_file_input_FreeFem,'w');
fprintf(fid,'%s\n',name_file_msh);
fprintf(fid,'%12.8f\n',lx);
fprintf(fid,'%12.8f\n',ly);
fprintf(fid,'%d\n',nx);
fprintf(fid,'%d\n',ny);
fclose(fid);


nb_layers=1;
multilayer(1).d=ly;
multilayer(1).mat=0;
% Termination condition // 0 for rigid backing 1 for radiation
termination=0;

% Number of waves (two ways) in each layer
% For the resolution, the incident waves is included in the system and put to RHS at the
% end of the procedure

% Addition of a new layer for the incident medium
l0.d=0;
l0.mat=0;
multilayer=[l0 multilayer];
nb_layers=nb_layers+1;

compute_number_PW_TMM
