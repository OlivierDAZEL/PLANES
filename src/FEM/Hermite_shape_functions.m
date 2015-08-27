function [Psi_1,Psi_2,Psi_3,Psi_4]=Hermite_shape_functions(lx)

Psi_1=[1 0/lx -3/lx/lx  2/lx/lx/lx];
Psi_2=[0 1/lx -2/lx/lx  1/lx/lx/lx];
Psi_3=[0 0/lx  -1/lx/lx 1/lx/lx/lx];
Psi_4=[0 0/lx  3/lx/lx -2/lx/lx/lx];

