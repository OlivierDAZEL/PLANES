% properties_PEM.m
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




% Biot densities with tortuosity effects

rho_12=-phi*air.rho*(alpha_tensor-eye(3));
rho_11=rho_1*eye(3)-rho_12;
rho_2=phi*air.rho*eye(3);
rho_22=rho_2-rho_12;

rho_22_til=phi^2*rho_eq_til;
rho_12_til=rho_2-rho_22_til;
rho_11_til=rho_1*eye(3)-rho_12_til;
rho_til=rho_11_til-((rho_12_til^2)*inv(rho_22_til));

gamma_til=phi*(rho_12_til*inv(rho_22_til)-((1-phi)/phi)*eye(3));
rho_s_til=rho_til+gamma_til^2*rho_eq_til;



switch porous_model.frame
    case{'anelastic'}
        structural_loss=1+(b_hat*(1i*omega/beta_hat)^alpha_hat)/(1+(1i*omega/beta_hat)^alpha_hat);
    case{'structural'}
        structural_loss=1+1j*eta;
end

if strcmp(porous_model.aniso,'yes')
    C_hat=C_hat_conservative*structural_loss;
    


  %  C_hat=rotate_S6(C_hat,Q');
        
%     C_hat=C_hat_conservative*structural_loss;
%     C_tot=C_tot_0;

%     rho_eq_til=Q*rho_eq_til*Q.';
%     rho_s_til= Q*rho_s_til *Q.';    
%     gamma_til= Q*gamma_til *Q.';

end