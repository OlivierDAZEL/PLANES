%   Program linear_triangle_acoustics.m
%
%   Copyright (C) 2014 LAUM UMR CNRS 6613 (France)
% 	   Olivier DAZEL <olivier.dazel@univ-lemans.fr>
%
% This program is free software: you can redistribute it and/or modify
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

%
% This function computes the H and Q elementary matrices of an acoustic linear element. 
%  It is a part of the course
%   "Numerical Methods in Acoustics and Vibration"
%  given for students of the MSc Master Program in Acoustics and 
%  IMDEA students
% 


function [vh,vq] = TR3_fluid(vcore)


c1 = 1/6;
a1=1/6;
b1=2/3;
ksig=[ a1, a1;
    a1, b1;
    b1,   a1];
cg=[c1,c1,c1];

npg=6;
a_gauss=0.445948490915965;
b_gauss=0.091576213509771;
c_gauss=0.111690794839005;
d_gauss=0.054975871827661;


ksig=[  a_gauss, a_gauss;
        1-2*a_gauss, a_gauss;
        a_gauss,   1-2*a_gauss;
        b_gauss, b_gauss;
        1-2*b_gauss, b_gauss;
        b_gauss, 1-2*b_gauss];
    
cg=[c_gauss,c_gauss,c_gauss,d_gauss,d_gauss,d_gauss];


% Matrix initialization
%-----------------------

vh=zeros(3,3);
vq=zeros(3,3);

% loop on Gauss Points
%------------------


for ipg=1:npg,

    ksi=ksig(ipg,1); 
    eta=ksig(ipg,2);
    N =  [ 1-ksi-eta  ksi  eta ];
    Nksi = [ -1 1 0 ];
    Neta = [ -1 0 1 ];
    dN = [ Nksi;Neta ];
    J = dN*vcore';
    j = inv(J);
    det = J(1,1)*J(2,2) - J(1,2)*J(2,1);
    poids = cg(ipg) * det;

    vB =j*dN;
   
    B_gdP=[vB(1,1) vB(1,2) vB(1,3);
           vB(2,1) vB(2,2) vB(2,3)];
    
    
    
    N_P=[N(1) N(2) N(3)];

    dh=B_gdP'*B_gdP*poids;
    vh=vh+dh;
    
    dq=N_P'*N_P*poids;
    vq=vq+dq;

    
end;   
