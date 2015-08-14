Psi_1_x=[1 0/lx -3/lx/lx  2/lx/lx/lx];
Psi_2_x=[0 1/lx -2/lx/lx  1/lx/lx/lx];
Psi_3_x=[0 0/lx  -1/lx/lx 1/lx/lx/lx];
Psi_4_x=[0 0/lx  3/lx/lx -2/lx/lx/lx];

Psi_1_y=[1 0/ly -3/ly/ly  2/ly/ly/ly];
Psi_2_y=[0 1/ly -2/ly/ly  1/ly/ly/ly];
Psi_3_y=[0 0/ly -1/ly/ly  1/ly/ly/ly];
Psi_4_y=[0 0/ly  3/ly/ly -2/ly/ly/ly];

Psi_1=multiply_polynom_2D(Psi_1_x,Psi_1_y');
Psi_2=multiply_polynom_2D(Psi_2_x,Psi_1_y');
Psi_3=multiply_polynom_2D(Psi_1_x,Psi_2_y');
Psi_4=multiply_polynom_2D(Psi_4_x,Psi_1_y');
Psi_5=multiply_polynom_2D(Psi_3_x,Psi_1_y');
Psi_6=multiply_polynom_2D(Psi_4_x,Psi_2_y');
Psi_7=multiply_polynom_2D(Psi_4_x,Psi_4_y');
Psi_8=multiply_polynom_2D(Psi_3_x,Psi_4_y');
Psi_9=multiply_polynom_2D(Psi_4_x,Psi_3_y');
Psi_10=multiply_polynom_2D(Psi_1_x,Psi_4_y');
Psi_11=multiply_polynom_2D(Psi_2_x,Psi_4_y');
Psi_12=multiply_polynom_2D(Psi_1_x,Psi_3_y');


