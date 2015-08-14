function [Psi_1,Psi_2,Psi_3,Psi_4]=Hermite_shape_functions_shifted(lx,x)

order_1=[-x 1]/lx;
order_2=multiply_polynom(order_1,order_1);
order_3=multiply_polynom(order_2,order_1);


Psi_1=add_polynom(add_polynom(1,-3*order_2),2*order_3);
Psi_2=add_polynom(add_polynom(order_1,-2*order_2),order_3);
Psi_3=add_polynom(-order_2,order_3);
Psi_4=add_polynom(3*order_2,-2*order_3);


