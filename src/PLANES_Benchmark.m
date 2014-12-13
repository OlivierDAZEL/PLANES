clear all
close all
clc
cd('Utils')
    init_path
cd ..

export_profiles=1;
nb_frequencies=1;


name_of_project='Benchmark1';
subproject=0;
plot_abs=1;
plot_TL=0;
init_PLANES

if plot_abs==1
    PLANES_result=load([name_project_directory name_of_project '.abs']);
    reference=load([name_project_directory name_of_project '.ref']);
    figure
    hold on
    plot(PLANES_result(:,1),PLANES_result(:,2),'r.')
    plot(reference(:,1),reference(:,4))
end
    
    