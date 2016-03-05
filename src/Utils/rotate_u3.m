function M=rotate_u3(angle)


% A       = cos(angle_x);
% B       = sin(angle_x);
% C       = cos(angle_y);
% D       = sin(angle_y);
% E       = cos(angle_z);
% F       = sin(angle_z);
% 
% AD      =   A * D;
% BD      =   B * D;
% 
% M(1,1)  =   C * E;
% M(1,2)  =  -C * F;
% M(1,3)  =  D;
% M(2,1)  =  BD * E + A * F;
% M(2,2)  =  -BD * F + A * E;
% M(2,3)  =  -B * C;
% M(3,1)  = -AD * E + B * F;
% M(3,2)  =  AD * F + B * E;
% M(3,3) =   A * C;


a_x=[1 0 0;0 cos(angle(1)) sin(angle(1));0 -sin(angle(1)) cos(angle(1))];
a_y=[cos(angle(2)) 0 -sin(angle(2));0 1 0;sin(angle(2)) 0 cos(angle(2))];
a_z=[cos(angle(3)) sin(angle(3)) 0;-sin(angle(3)) cos(angle(3)) 0;0 0 1];

M=a_x*a_y*a_z;





















