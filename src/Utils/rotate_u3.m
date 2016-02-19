function M=rotate_u3(angle_x,angle_y,angle_z)


A       = cos(angle_x);
B       = sin(angle_x);
C       = cos(angle_y);
D       = sin(angle_y);
E       = cos(angle_z);
F       = sin(angle_z);

AD      =   A * D;
BD      =   B * D;

M(1,1)  =   C * E;
M(1,2)  =  -C * F;
M(1,3)  =  D;
M(2,1)  =  BD * E + A * F;
M(2,2)  =  -BD * F + A * E;
M(2,3)  =  -B * C;
M(3,1)  = -AD * E + B * F;
M(3,2)  =  AD * F + B * E;
M(3,3) =   A * C;

