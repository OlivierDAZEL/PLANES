clear all
close all
current_dir=pwd;


frequency.nb=-100;
frequency.min=1;
frequency.max=20000;


multilayer(1,1).nb=1;
multilayer(1,1).d=5e-2;
multilayer(1,1).mat=5010;
% multilayer(2,1).d=1e-2;
% multilayer(2,1).mat=5003;
% multilayer(3,1).d=1e-2;
% multilayer(3,1).mat=5003;
% Termination condition // 0 for rigid backing 1 for radiation
multilayer(1,1).termination=0;

% multilayer(1,2).nb=2;
% multilayer(1,2).d=1e-2;
% multilayer(1,2).mat=5003;
% multilayer(2,2).d=1e-2;
% multilayer(2,2).mat=5003;
% % Termination condition // 0 for rigid backing 1 for radiation
% multilayer(1,2).termination=0;

data_model.theta_inc=45*pi/180;
data_model.theta_x=23*pi/180;
data_model.theta_y=0*pi/180;
cd ../../src/Main/
PLANES_Multilayer_3D%('Multilayer_3D',101,data_model,multilayer,frequency)

cd(current_dir);

load('test.mat');
test1_JPPM=load('out/Multilayer_3D_101.PW3D');



figure 
semilogx(test(1).f,test(1).absorp)
hold on
semilogx(test1_JPPM(:,1),test1_JPPM(:,2),'r')
