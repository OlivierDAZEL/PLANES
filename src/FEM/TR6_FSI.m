% TR6_FSI.m
%
% Copyright (C) 2014 < Olivier DAZEL <olivier.dazel@univ-lemans.fr> >
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


function FSIe = TR6_FSI(a1,a2)    J=norm(a2-a1)/2.00D+00;            npg=4;    ksi_g(1)= 0.339981043584856;    ksi_g(2)=-0.339981043584856;    ksi_g(3)= 0.861136311594053;    ksi_g(4)=-0.861136311594053;          p_g(1)=   0.652145154862546;    p_g(2)=   0.652145154862546;    p_g(3)=   0.347854845137454;    p_g(4)=   0.347854845137454;        FSIe=0.00D+00;	        for ipg=1:npg        	ksi=ksi_g(ipg);    	    	N(1,1)=-(1-ksi)*ksi/2;    	N(2,1)= (1+ksi)*ksi/2;    	N(3,1)= (1-ksi)*(1+ksi);    	    	FSIe=FSIe+N*N'*J*p_g(ipg);     end    