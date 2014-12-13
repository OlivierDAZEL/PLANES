warning('off','MATLAB:nearlySingularMatrix')


name_project_directory=['../Projects/' name_of_project '/'];

if subproject==0
    name_file=[name_project_directory name_of_project];
    name_of_project_full=name_of_project;
else
    name_file=[name_project_directory name_of_project '_' num2str(subproject) ];
    name_of_project_full=[name_of_project '_' num2str(subproject)];
end

name_file_msh=          [name_file '.msh'];
name_file_edp=          [name_file '.edp'];
name_file_input_FreeFem=['FF.inp'];
name_file_abs=          [name_file '.abs'];
name_file_TL=           [name_file '.TL'];
name_file_info=          [name_file '.info'];

if export_profiles==1
    name_directory_profiles= [name_project_directory '/Profiles/'];
    if ~exist(name_directory_profiles,'dir')
        mkdir(name_directory_profiles)
    end
end



set(0,'DefaultLineMarkerSize',15);
set(0,'Defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultAxeslinewidth',2);



