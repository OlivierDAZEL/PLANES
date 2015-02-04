function M=Phi_fluid_2D(nx,ny,Z_0)

M(1:3,1:length(nx))=[-nx;-ny;Z_0*ones(1,length(nx))];