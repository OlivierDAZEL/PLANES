clear all
close all
current_dir=pwd;
PLANES_dir='../../src/Main/';
cd(PLANES_dir)


frequency.nb=1;
frequency.min=200;
frequency.max=20000;

 PLANES_Multilayer('Multilayer_3D',107,data_model,multilayer_6,frequency)

  cd(current_dir);

