function [Psi_1,Psi_2,Psi_3,Psi_4,Psi_5,Psi_6,Psi_7,Psi_8,Psi_9,Psi_10,Psi_11,Psi_12]=H12_shape_functions(lx,ly)

[Psi_1_x,Psi_2_x,Psi_3_x,Psi_4_x]=Hermite_shape_functions(lx);
[Psi_1_y,Psi_2_y,Psi_3_y,Psi_4_y]=Hermite_shape_functions(ly);

Psi_1 =multiply_polynom_2D(Psi_1_x,Psi_1_y');
Psi_2 =multiply_polynom_2D(Psi_2_x,Psi_1_y');
Psi_3 =multiply_polynom_2D(Psi_1_x,Psi_2_y');
Psi_4 =multiply_polynom_2D(Psi_4_x,Psi_1_y');
Psi_5 =multiply_polynom_2D(Psi_3_x,Psi_1_y');
Psi_6 =multiply_polynom_2D(Psi_4_x,Psi_2_y');
Psi_7 =multiply_polynom_2D(Psi_4_x,Psi_4_y');
Psi_8 =multiply_polynom_2D(Psi_3_x,Psi_4_y');
Psi_9 =multiply_polynom_2D(Psi_4_x,Psi_3_y');
Psi_10=multiply_polynom_2D(Psi_1_x,Psi_4_y');
Psi_11=multiply_polynom_2D(Psi_2_x,Psi_4_y');
Psi_12=multiply_polynom_2D(Psi_1_x,Psi_3_y');


