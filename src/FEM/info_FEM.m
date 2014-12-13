fidinfo=fopen(name_file_info,'w');
fprintf(fidinfo,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
fprintf(fidinfo,'!!                  Output files of PLANES                   !!\n');
fprintf(fidinfo,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
[~,name]=system('hostname')
fprintf(fidinfo,'!! Generated %s on %s', datestr(now,'dd-mm-yyyy at HH-MM-SS'),name);
fprintf(fidinfo,'!! Name of project = %s\n',name_of_project);
fprintf(fidinfo,'!! Subproject # %d\n',subproject);
fprintf(fidinfo,'!! #dof = %d\n',length(X));
fprintf(fidinfo,'!! Computation time = %d s\n',time_FEM);
fclose(fidinfo)