clear all
close all
current_dir=pwd;


frequency.nb=-200;
frequency.min=1;
frequency.max=20000;


multilayer(1,1).nb=1;
multilayer(1,1).d=5e-2;
multilayer(1,1).mat=5010;
% Termination condition // 0 for rigid backing 1 for radiation
multilayer(1,1).termination=0;

% multilayer(1,2).nb=2;
% multilayer(1,2).d=1e-2;
% multilayer(1,2).mat=5010;
% multilayer(2,2).d=4e-2;
% multilayer(2,2).mat=5010;
% multilayer(1,2).termination=0;
% 
% multilayer(1,3).nb=3;
% multilayer(1,3).d=1e-2;
% multilayer(1,3).mat=5010;
% multilayer(2,3).d=3e-2;
% multilayer(2,3).mat=5010;
% multilayer(3,3).d=1e-2;
% multilayer(3,3).mat=5010;
% multilayer(1,3).termination=0;

data_model.theta=[23 62]*pi/180;

cd ../../src/Main/
PLANES_Multilayer('Multilayer_3D',101,data_model,multilayer,frequency)

cd(current_dir);

load('160303_test8.mat')
test1_JPPM=load('out/Multilayer_3D_101.PW');


figure 
semilogx(test1_JPPM(:,1),absorp)
hold on
semilogx(test1_JPPM(:,1),test1_JPPM(:,2),'r.')
% semilogx(test1_JPPM(:,1),test1_JPPM(:,6),'m+')
% semilogx(test1_JPPM(:,1),test1_JPPM(:,10),'k')
