clear all
close all

l_x=2
x_shift=0.5
nb_points=20;

[Psi_1,Psi_2,Psi_3,Psi_4]=Hermite_shape_functions(l_x);
[Psi_1_shift,Psi_2_shift,Psi_3_shift,Psi_4_shift]=Hermite_shape_functions_shifted(l_x,x_shift);


x=linspace(0,l_x,nb_points);

for ii=1:nb_points
   x_temp=x(ii); 
   H1(ii)=evaluate_polynom_1D(Psi_1,x_temp); 
   H4(ii)=evaluate_polynom_1D(Psi_4,x_temp);
   x_temp=x(ii)+x_shift;
   H1_shift(ii)=evaluate_polynom_1D(Psi_1_shift,x_temp); 
   H4_shift(ii)=evaluate_polynom_1D(Psi_4_shift,x_temp);   
   
end

figure
hold on
plot(x,H1,'b')
plot(x,H4,'r')
plot(x,H1_shift,'b.')
plot(x,H4_shift,'r.')