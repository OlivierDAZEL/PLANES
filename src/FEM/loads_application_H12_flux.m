% loads_application_H12_flux.m
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

for ie=1:nb.loads
    typ=floor(loads(ie,4));
    
    
    
    node_1=loads(ie,1);
    node_2=loads(ie,2);
    
    x1=nodes(loads(ie,1),1);
    y1=nodes(loads(ie,1),2);
    x2=nodes(loads(ie,2),1);
    y2=nodes(loads(ie,2),2);

    
    
    length_x=sqrt((x2-x1)^2+(y2-y1)^2);

    load_Hermite_2D
    index_p=dof_A(p_H12(elements(loads(ie,3),:)));
    indice_test=index_p([1 2 5 4]);
    indice_champs_0=indice_test;
    indice_champs_y=index_p([3 6]);
   
    
    
    switch typ
        case {3}
            if sort(loads(ie,1:2))==sort(elements(loads(ie,3),1:2)) % Interface node 1 node 2: bottom edge
                C=[0 1 0];
                s=1;
                
                
                
                
                
                
                
                
                
                nx=0;
                ny=-1;
                
                

                W=[nx -ny -nx;ny nx -ny;air.Z 0 air.Z];
                Omega=[nx/2 ny/2 1/(2*air.Z);-ny nx 0;-nx/2 -ny/2 1/(2*air.Z)];
                
                W0plus=W(:,1:2);
                Wmoins=W(:,3);
                Omega0plus=Omega(1:2,:);
                Omegamoins=Omega(3,:);
                S_tilde=inv(C*Wmoins)*s;
                R_tilde=inv(C*Wmoins)*(-C*W0plus);
                FF=Wmoins*S_tilde;
                PP=(W0plus+Wmoins*R_tilde)*Omega0plus;
                
                line_PP=2;
                PP=-PP(line_PP,:)/(1j*omega); %\d p/\p n =-v_y
                FF=-FF(line_PP)/(1j*omega);
                for i_test=1:4
                    eval(['Psi_test=Psi_',num2str(i_test),'_x;'])
                    F(indice_tesedit(i_test))=F(indice_test(i_test))+integrate_polynom(Psi_test,lx_H12)*FF;
                    for i_champs=1:4
                        eval(['temp_1=PP(3)*Psi_',num2str(i_champs),'_x;']);
                        eval(['temp_2=PP(1)*(-1/(j*air.rho*omega))*derive_polynom(Psi_',num2str(i_champs),'_x);']);
                        temp=add_polynom(temp_1,temp_2);
                        A(indice_test(i_test),indice_champs_0(i_champs))=A(indice_test(i_test),indice_champs_0(i_champs))-integrate_polynom(multiply_polynom(Psi_test,temp),lx_H12);
                    end
                    for i_champs=1:2
                        index_Psi=[1 4];
                        %['temp=(1/ly_H12)*PP(2)*(-1/(j*air.rho*omega))*Psi_',num2str(index_Psi(i_champs)),'_x;']
                        eval(['temp=(1/ly_H12)*PP(2)*(-1/(j*air.rho*omega))*Psi_',num2str(index_Psi(i_champs)),'_x;'])
                        A(indice_test(i_test),indice_champs_y(i_champs))=A(indice_test(i_test),indice_champs_y(i_champs))-integrate_polynom(multiply_polynom(Psi_test,temp),lx_H12);
                    end
                end
            else 
                aezezaezaezeazezaeazezaeazezaeza
            end
            
    end
end


