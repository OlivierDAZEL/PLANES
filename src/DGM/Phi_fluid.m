function M=Phi_fluid_2D(nx,ny,Z_e)

M(1:3,1:length(nx))=[-nx;-ny;Z_e*ones(1,length(nx))];