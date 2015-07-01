Psi_1_x=[1 0/lx_H12 -3/lx_H12/lx_H12  2/lx_H12/lx_H12/lx_H12];
Psi_2_x=[0 1/lx_H12 -2/lx_H12/lx_H12  1/lx_H12/lx_H12/lx_H12];
Psi_3_x=[0 0/lx_H12  -1/lx_H12/lx_H12 1/lx_H12/lx_H12/lx_H12];
Psi_4_x=[0 0/lx_H12  3/lx_H12/lx_H12 -2/lx_H12/lx_H12/lx_H12];

Psi_1_y=[1 0/ly_H12 -3/ly_H12/ly_H12  2/ly_H12/ly_H12/ly_H12];
Psi_2_y=[0 1/ly_H12 -2/ly_H12/ly_H12  1/ly_H12/ly_H12/ly_H12];
Psi_3_y=[0 0/ly_H12 -1/ly_H12/ly_H12  1/ly_H12/ly_H12/ly_H12];
Psi_4_y=[0 0/ly_H12  3/ly_H12/ly_H12 -2/ly_H12/ly_H12/ly_H12];

Psi_11=multiply_polynom_2D(Psi_1_x,Psi_1_y');
Psi_12=multiply_polynom_2D(Psi_1_x,Psi_2_y');
Psi_13=multiply_polynom_2D(Psi_1_x,Psi_3_y');
Psi_14=multiply_polynom_2D(Psi_1_x,Psi_4_y');
Psi_21=multiply_polynom_2D(Psi_2_x,Psi_1_y');
Psi_22=multiply_polynom_2D(Psi_2_x,Psi_2_y');
Psi_23=multiply_polynom_2D(Psi_2_x,Psi_3_y');
Psi_24=multiply_polynom_2D(Psi_2_x,Psi_4_y');
Psi_31=multiply_polynom_2D(Psi_3_x,Psi_1_y');
Psi_32=multiply_polynom_2D(Psi_3_x,Psi_2_y');
Psi_33=multiply_polynom_2D(Psi_3_x,Psi_3_y');
Psi_34=multiply_polynom_2D(Psi_3_x,Psi_4_y');
Psi_41=multiply_polynom_2D(Psi_4_x,Psi_1_y');
Psi_42=multiply_polynom_2D(Psi_4_x,Psi_2_y');
Psi_43=multiply_polynom_2D(Psi_4_x,Psi_3_y');
Psi_44=multiply_polynom_2D(Psi_4_x,Psi_4_y');

Xi(1,1,:,:)=Psi_11
qsdsqdsqsqd