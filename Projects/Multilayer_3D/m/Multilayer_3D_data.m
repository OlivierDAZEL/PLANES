
if (nargin==0)
    
    frequency.nb=100;
    frequency.min=50;
    frequency.max=5000;
    
    profiles.mesh=0;
    profiles.solution=0;
    profiles.x=0;
    profiles.y=0;
    profiles.map=0;
    profiles.custom=0;
    profiles.on=profiles.x+profiles.y+profiles.map+profiles.custom;
    profiles.custom_plots = {};
    export.profiles=0;
    export.nrj=1;
    
    data_model.thicknessplate1=1e-3;
    data_model.thicknessporous=2e-2;
    data_model.thicknessplate2=1e-3;
    
    data_model.labelplate1=1001;
    data_model.labelplate2=1001;
    data_model.labelporous=5003;
    
    data_model.theta_inc=45*pi/180;
    data_model.theta_x=data_model.theta_inc;
    data_model.theta_y=0;
    
    
end


nb_multilayers=2;

nb_layers(1)=2;
multilayer(1,1).d=data_model.thicknessporous/2;
multilayer(1,1).mat=data_model.labelporous;
multilayer(2,1).d=data_model.thicknessporous/2;
multilayer(2,1).mat=data_model.labelporous;
% Termination condition // 0 for rigid backing 1 for radiation
termination(1)=0;

nb_layers(2)=1;
multilayer(1,2).d=data_model.thicknessporous;
multilayer(1,2).mat=data_model.labelporous;
% Termination condition // 0 for rigid backing 1 for radiation
termination(2)=0;


nb_multilayers_3D=nb_multilayers;
nb_layers_3D=nb_layers;
multilayer_3D=multilayer;
termination_3D=termination;

%clear multilayer

