clear all
close all
current_dir=pwd;
PLANES_dir='../../src/Main/';
cd(PLANES_dir)


frequency.nb=-200;
frequency.min=1;
frequency.max=20000;

% Test 1
multilayer_1(1,1).nb=1;
multilayer_1(1,1).d=5e-2;
multilayer_1(1,1).mat=5005;
multilayer_1(1,1).termination=0;

data_model.theta_inc=[23]*pi/180
% data_model.theta_inc=[23 62]*pi/180;
% 
 PLANES_PW('Multilayer_3D',101,data_model,multilayer_1,frequency)
% 
% 
% 
% % Test 2
% multilayer_2=multilayer_1;
% multilayer_2(1,1).termination=1;
% PLANES_Multilayer('Multilayer_3D',102,data_model,multilayer_2,frequency)
% % 
% % % Test 3
% % data_model.theta=[1e-3 1e-3]*pi/180;
% % multilayer_3(1,1).nb=3;
% % multilayer_3(1,1).d=1e-3;
% % multilayer_3(1,1).mat=1016;
% % multilayer_3(2,1).d=88e-3;
% % multilayer_3(2,1).mat=5012;
% % multilayer_3(3,1).d=1e-3;
% % multilayer_3(3,1).mat=1016;
% % multilayer_3(1,1).termination=1;
% % PLANES_Multilayer('Multilayer_3D',103,data_model,multilayer_3,frequency)
% % % 
% % % Test 4
% % 
% % multilayer_4=multilayer_3;
% % data_model.theta=[45 50]*pi/180;
% % PLANES_Multilayer('Multilayer_3D',104,data_model,multilayer_4,frequency)
% % 
% % % Test 5
% % multilayer_5=multilayer_4;
% % multilayer_5(2,1).mat=5013;
% % PLANES_Multilayer('Multilayer_3D',105,data_model,multilayer_5,frequency)
% % 
% % % 
% % % multilayer_6(1,1).nb=4;
% % % multilayer_6(1,1).d=1e-3;
% % % multilayer_6(1,1).mat=1016;
% % % multilayer_6(2,1).d=88e-3;
% % % multilayer_6(2,1).mat=5012;
% % % multilayer_6(3,1).d=88e-3;
% % % multilayer_6(3,1).mat=5013;
% % % multilayer_6(4,1).d=1e-3;
% % % multilayer_6(4,1).mat=1016;
% % % multilayer_6(1,1).termination=1;
% % % 
% % % PLANES_Multilayer('Multilayer_3D',106,data_model,multilayer_6,frequency)
% % % 
% % % multilayer_6(2,1).mat=5013;
% % % multilayer_6(3,1).mat=5012;
% % % 
% % % PLANES_Multilayer('Multilayer_3D',107,data_model,multilayer_6,frequency)
% % % 
% % % 
%   cd(current_dir);
% % 
%  load('160311_tests.mat')
% % 
%  test1_JPPM=load('out/Multilayer_3D_101.PW');
%  test2_JPPM=load('out/Multilayer_3D_102.PW');
% % test3_JPPM=load('out/Multilayer_3D_103.PW');
% %test4_JPPM=load('out/Multilayer_3D_104.PW');
% % test5_JPPM=load('out/Multilayer_3D_105.PW');
% % % test6_JPPM=load('out/Multilayer_3D_106.PW');
% % % test7_JPPM=load('out/Multilayer_3D_107.PW');
% % 
% figure 
% semilogx(test(1).f,test(1).absorp)
% hold on
% semilogx(test1_JPPM(:,1),test1_JPPM(:,2),'r.')
% % 
% figure 
% semilogx(test(2).f,test(2).TL)
% hold on
% semilogx(test2_JPPM(:,1),test2_JPPM(:,5),'r.')
% % 
% figure 
% semilogx(test(3).f,test(3).TL)
% hold on
% semilogx(test3_JPPM(:,1),test3_JPPM(:,5),'r.')
% % 
% % 
% figure 
% semilogx(test(4).f,test(4).TL)
% hold on
% semilogx(test4_JPPM(:,1),test4_JPPM(:,5),'r.')
% % 
% % 
% % figure 
% % semilogx(test(5).f,test(5).TL)
% % hold on
% % semilogx(test5_JPPM(:,1),test5_JPPM(:,5),'r.')
% % 
% % % figure 
% % % semilogx(test6_JPPM(:,1),TL6,'b')
% % % hold on
% % % semilogx(test6_JPPM(:,1),TL7,'r')
% % % semilogx(test6_JPPM(:,1),test6_JPPM(:,5),'b.')
% % % semilogx(test7_JPPM(:,1),test7_JPPM(:,5),'r.')
% % % 
% % % figure 
% % % semilogx(test(7).f,test(7).TL)
% % % hold on
% % % semilogx(test7_JPPM(:,1),test7_JPPM(:,5),'r.')
% % % 
