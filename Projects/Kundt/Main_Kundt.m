clear all
close all
current_dir=pwd;
PLANES_dir='../../src/Main/';
cd(PLANES_dir)


frequency.nb=1;
frequency.min=500;
frequency.min=100;
frequency.max=20000;

data_model.profiles.mesh=1;
data_model.profiles.solution=0;
data_model.profiles.x=0;
data_model.profiles.y=1;
data_model.profiles.map=0;
data_model.profiles.custom=0;
%data_model.profiles.on=data_model.profiles.x+data_model.profiles.y+data_model.profiles.map;
data_model.profiles.custom_plots = {};
data_model.export.profiles=0;
data_model.export.nrj=1;

data_model.lx=0.1;
data_model.ly=1;
data_model.nx=1;
%data_model.ny=ceil(data_model.nx*data_model.ly/data_model.lx);
data_model.ny=20;

data_model.theta_DGM.nb=4;
data_model.tilt=0*pi/data_model.theta_DGM.nb;


PLANES('Kundt',7,data_model,frequency)

cd(current_dir);

% Computation of the analytical solution and superposition

cd('../../src/Physics')
air_properties_maine
cd(current_dir);
omega=2*pi*frequency.min;
k_air=omega/air.c;

% u=A*sin(k_air(x-data_model.ly))
% 1=A*sin(-k_air*data_model.ly)
A=1/(sin(-k_air*data_model.ly)*(1j*omega));

x=linspace(0,data_model.ly,200);

p=-air.K*k_air*A*cos(k_air*(x-data_model.ly));
figure(2010)
hold on
plot(x,abs(p),'b')
figure(4010)
hold on
plot(x,angle(p),'b')
