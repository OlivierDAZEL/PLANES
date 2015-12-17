project.name='Example';
project.num=1;

% Number of frequencies
% If the number is negative then a logscale is chosen
% If the number is equal to 1, then the frequency is equal to freq_min
frequency.nb=1;
frequency.min=100;
frequency.max=750;

profiles.mesh=1;
profiles.solution=0;
profiles.x=0;
profiles.y=0;
profiles.map=1;
profiles.custom=0;
profiles.on=profiles.x+profiles.y+profiles.map+profiles.custom;
profiles.custom_plots = {};

export.nrj=0;
export.profiles=0;
