period=55e-3;
thicknessporous=5.5e-2;
labelporous=5005;
labelinclusion=1004;
radius=1.25e-2;
thicknessinclusion=2e-3;
nx=20;
nepaisseurinclusion=1;
ninclusion=ceil(nepaisseurinclusion*2*pi*radius/thicknessinclusion);


fid=fopen(name_file_input_FreeFem,'w');
fprintf(fid,'%s\n',name_file_msh);
fprintf(fid,'%12.8f\n',period);
fprintf(fid,'%12.8f\n',thicknessporous);
fprintf(fid,'%d\n',labelporous);
fprintf(fid,'%d\n',labelinclusion);
fprintf(fid,'%12.8f\n',radius);
fprintf(fid,'%12.8f\n',thicknessinclusion);
fprintf(fid,'%d\n',nx);
fprintf(fid,'%d\n',nepaisseurinclusion);
fprintf(fid,'%d\n',ninclusion);
fclose(fid);


nb_layers=1;
multilayer(1).d=thicknessporous;
multilayer(1).mat=labelporous;
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
