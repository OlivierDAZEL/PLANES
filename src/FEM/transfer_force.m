% transfer_force.m
%
% Copyright (C) 2015 < Olivier DAZEL <olivier.dazel@univ-lemans.fr> >
%
% This file is part of PLANES.
%
% PLANES (Porous LAum NumErical Simulator) is a software to compute the
% vibroacoustic response of sound packages containing coupled
% acoustic/elastic/porous substructures. It is mainly based on the
% Finite-Element Method and some numerical methods developped at
% LAUM (http://laum.univ-lemans.fr).
%
% You can download the latest version of PLANES at
% https://github.com/OlivierDAZEL/PLANES
% or find more details on Olivier's webpage
% http://perso.univ-lemans.fr/~odazel/
%
% For any questions or if you want to
% contribute to this project, contact Olivier.
%
% PLANES is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%%


function    F_plus=transfer_force(k_x,omega,F_moins,media)

lambda=media(1).lambda;
mu=media(1).mu;
rho=media(1).rho;
d=media(1).thickness;

P=lambda+2*mu;

delta_P=omega*sqrt(rho/P);
delta_s=omega*sqrt(rho/mu);


beta_P=sqrt(delta_P^2-k_x^2);
beta_s=sqrt(delta_s^2-k_x^2);


alpha_P=-1j*lambda*delta_P^2-1j*2*mu*beta_P^2;
alpha_s= 2*1j*mu*beta_s*k_x;

V_0=diag([1j*beta_P,-1j*beta_P,1j*beta_s,-1j*beta_s]);


Phi_0(1,1)=-2*1j*mu*beta_P*k_x;
Phi_0(1,2)=2*1j*mu*beta_P*k_x;
Phi_0(1,3)=1j*mu*(beta_s^2-k_x^2);
Phi_0(1,4)=1j*mu*(beta_s^2-k_x^2);


Phi_0(2,1)= beta_P;
Phi_0(2,2)=-beta_P;
Phi_0(2,3)= k_x;
Phi_0(2,4)= k_x;

Phi_0(3,1)=alpha_P;
Phi_0(3,2)=alpha_P;
Phi_0(3,3)=-alpha_s;
Phi_0(3,4)=alpha_s;

Phi_0(4,1)=k_x;
Phi_0(4,2)=k_x;
Phi_0(4,3)=-beta_s;
Phi_0(4,4)=beta_s;

[a,indice]=sort(diag(real(V_0)),'descend');

for i_m=1:4
    Phi(:,i_m)=Phi_0(:,indice(4+1-i_m));
    lamda(i_m)=-V_0(indice(4+1-i_m),indice(4+1-i_m));
end

Psi=inv(Phi);

if k_x~=0
    
    
    MM=exp(lamda(1)*d)*(Phi(:,1)*Psi(1,:)+exp((lamda(2)-lamda(1))*d)*Phi(:,2)*Psi(2,:)+exp((lamda(3)-lamda(1))*d)*Phi(:,3)*Psi(3,:)+exp((lamda(4)-lamda(1))*d)*Phi(:,4)*Psi(4,:));
    
    
    X_plus=exp(lamda(1)*d)*Psi(1,:)*F_moins;
    
    F_plus=Phi(:,1)+(exp((lamda(2)-lamda(1))*d)*Phi(:,2)*Psi(2,:)+exp((lamda(3)-lamda(1))*d)*Phi(:,3)*Psi(3,:)+exp((lamda(4)-lamda(1))*d)*Phi(:,4)*Psi(4,:))*F_moins/(Psi(1,:)*F_moins);
    
    F_plus=F_plus*X_plus;
    
    
else  % k_x==0
    
    MM=(Phi(:,1)*Psi(1,:)*exp(lamda(1)*d)+exp((lamda(2))*d)*Phi(:,2)*Psi(2,:)+exp((lamda(3))*d)*Phi(:,3)*Psi(3,:)+exp((lamda(4))*d)*Phi(:,4)*Psi(4,:));
    F_plus=MM*F_moins;
    
end