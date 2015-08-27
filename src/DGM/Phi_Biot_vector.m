% Phi_Biot_vector.m
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
%                 https://github.com/OlivierDAZEL/PLANES
% or find more details on Olivier's webpage
%                  http://perso.univ-lemans.fr/~odazel/
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
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
function M=Phi_Biot_vector(nx,ny,delta_1,delta_2,delta_3,mu_1,mu_2,mu_3,N,A_hat,K_eq_til,omega,Shift)

n=length(nx);

M=[nx;ny ;mu_1*nx;mu_1*ny ;-(delta_1/omega)*(A_hat+2*N*nx.^2);   -(delta_1/omega)*(2*N*nx.*ny)    ; -(delta_1/omega)*(A_hat+2*N*ny.^2);-(delta_1/omega)*(K_eq_til*mu_1)*ones(1,length(nx))];
M=[M [nx;ny ;mu_2*nx;mu_2*ny ;-(delta_2/omega)*(A_hat+2*N*nx.^2);-(delta_2/omega)*(2*N*nx.*ny)    ; -(delta_2/omega)*(A_hat+2*N*ny.^2);-(delta_2/omega)*(K_eq_til*mu_2)*ones(1,length(nx))]];
M=[M [ny;-nx;mu_3*ny;-mu_3*nx;-(delta_3/omega)*(2*N*nx.*ny)     ; (delta_3/omega)*(N*(nx.^2-ny.^2)); (delta_3/omega)*(2*N*nx.*ny)                         ;0*ones(1,length(nx))]];


% M(3,4)=ny;
% M(4,4)=-nx;
%
% M(5,5)=ny^2;
% M(6,5)=-nx*ny;
% M(7,5)=nx^2;

M5=(M(5,:)+M(7,:))/2;
M7=(M(5,:)-M(7,:))/2;

M(5,:)=M5;
M(7,:)=M7;

M=repmat(M,3*n,1);

M=Shift.*M;

