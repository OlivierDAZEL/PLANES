function [A_x_out,A_y_out,A_z_out]=rotate_Axyz(A_x,A_y,A_z,Q)

A_x_out=A_x;
A_y_out=A_y;
A_z_out=A_z;


rot=Q;

A_x_out(1:3,:)=rot(1,1)*A_x(1:3,:)+rot(1,2)*A_y(1:3,:)+rot(1,3)*A_z(1:3,:);
A_y_out(1:3,:)=rot(2,1)*A_x(1:3,:)+rot(2,2)*A_y(1:3,:)+rot(2,3)*A_z(1:3,:);
A_z_out(1:3,:)=rot(3,1)*A_x(1:3,:)+rot(3,2)*A_y(1:3,:)+rot(3,3)*A_z(1:3,:);


A_x_out(4:6,:)=rot(1,1)*A_x(4:6,:)+rot(1,2)*A_y(4:6,:)+rot(1,3)*A_z(4:6,:);
A_y_out(4:6,:)=rot(2,1)*A_x(4:6,:)+rot(2,2)*A_y(4:6,:)+rot(2,3)*A_z(4:6,:);
A_z_out(4:6,:)=rot(3,1)*A_x(4:6,:)+rot(3,2)*A_y(4:6,:)+rot(3,3)*A_z(4:6,:);



