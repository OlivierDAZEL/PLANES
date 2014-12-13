
thicknessplate=1e-3;
thicknessporous=1e-3;
labelplate=1001;
labelporous=5003;

fid=fopen(name_file_input_FreeFem,'w');
fprintf(fid,'%s\n',name_file_msh);
fclose(fid);


nb_layers=3;
multilayer(1).d=thicknessplate;
multilayer(1).mat=labelplate;
multilayer(2).d=thicknessporous;
multilayer(2).mat=labelporous;
multilayer(3).d=labelplate;
multilayer(3).mat=thicknessplate;
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