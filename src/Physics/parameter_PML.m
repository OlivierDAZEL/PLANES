function [tau_x,tau_y]=parameter_PML(ielem)


temp=(ielem-8000);
if (floor(temp/100)==1)
    tau_x=exp(1j*pi/4);
else
    tau_x=1;
end
temp=temp-100*floor(temp/100);
if (floor(temp/10)==1)
    tau_y=exp(1j*pi/4);
else
    tau_y=1;
end