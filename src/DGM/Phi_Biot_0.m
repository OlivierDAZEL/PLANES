function M=Phi_Biot_0(nx,ny)


M=zeros(8,2);

M(3,1)=ny;
M(4,1)=-nx;

M(5,2)=ny^2;
M(6,2)=-nx*ny;
M(7,2)=nx^2;

M5=(M(5,:)+M(7,:))/2;
M7=(M(5,:)-M(7,:))/2;

M(5,:)=M5;
M(7,:)=M7;