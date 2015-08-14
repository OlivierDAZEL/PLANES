% loads_application_TR6.m
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

for ie=1:nb.loads
    typ=floor(edges.loads(ie,4));
    
    x1=nodes(edges.loads(ie,1),1);
    y1=nodes(edges.loads(ie,1),2);
    x2=nodes(edges.loads(ie,2),1);
    y2=nodes(edges.loads(ie,2),2);
    length_edge=sqrt((x2-x1)^2+(y2-y1)^2);
    
    switch typ
        case {3}  % Unit velocity (classic)
            
            if elem.model(edges.loads(ie,3))==1
                if (x1<x2)
                    a=x1;
                    node(1)=edges.loads(ie,1);
                    node(2)=edges.loads(ie,2);
                    node(3)=edges.loads(ie,6);
                else
                    a=x2;
                    node(2)=edges.loads(ie,1);
                    node(1)=edges.loads(ie,2);
                    node(3)=edges.loads(ie,6);
                end
                F3=TR6_unit(length_edge);
                index_force=dof_A(p_TR(node));
                index_F_elem=find(index_force);
                index_F_global=index_force(index_F_elem);
                F(index_F_global)=F(index_F_global)-F3(index_F_elem)/(1j*omega);
            elseif  elem.model(edges.loads(ie,3))==2
                lx=norm(nodes(elem.nodes(edges.loads(ie,3),1),:)-nodes(elem.nodes(edges.loads(ie,3),2),:));
                ly=norm(nodes(elem.nodes(edges.loads(ie,3),1),:)-nodes(elem.nodes(edges.loads(ie,3),4),:));
                load_Hermite_2D_2
                index_p=dof_A(p_H12(elem.nodes(edges.loads(ie,3),1:4)));
                % What is the edge on the element ?
                if sort(edges.loads(ie,1:2))==sort(elem.nodes(edges.loads(ie,3),1:2)) % Bottom
                    indice_test=index_p([1 2 5 4]);
                elseif sort(edges.loads(ie,1:2))==sort(elem.nodes(edges.loads(ie,3),2:3))
                    indice_test=index_p([4 6 9 7]);
                elseif sort(edges.loads(ie,1:2))==sort(elem.nodes(edges.loads(ie,3),2:3))
                    indice_test=index_p([10 11 8 7]);
                elseif sort(edges.loads(ie,1:2))==sort(elem.nodes(edges.loads(ie,3),2:3))
                    indice_test=index_p([1 3 12 10]);
                end
                for i_test=1:4
                    eval(['Psi_test=Psi_',num2str(i_test),'_x;'])
                    F(indice_test(i_test))=F(indice_test(i_test))+integrate_polynom(Psi_test,length_edge)/(1j*omega);
                end
            elseif  ismember(elem.model(edges.loads(ie,3)),[10 11])
                boundary_normal_displacement
                
            end
        case {60}
            if sort(edges.loads(ie,1:2))==sort(elem.nodes(edges.loads(ie,3),1:2)) % Interface node 1 node 2: bottom edge
                
                
                c_e=air.c;
                k_e=omega/c_e;
                Z_e=air.Z;
                rho_e=air.rho;
                M_e=diag([air.rho,air.rho,1/air.K]);
                
                e_1=edges.loads(ie,3);
                c_1=mean(nodes(nonzeros(elem.nodes(e_1,:)),1:2))';
                e_edge=e_1;
                
                
                coord_edge(1:2,1)=nodes(edges.loads(ie,1),1:2)';
                coord_edge(1:2,2)=nodes(edges.loads(ie,2),1:2)';
                
                a=coord_edge(:,1);
                b=coord_edge(:,2);
                
                h=norm(b-a);
                n_ab=(b-a)/h;
                
                
                %%%%% vector normal pointing out from e
                
                centre_edge=(a+b)/2;
                n_centre=c_1-centre_edge;
                ne=normal_edge(coord_edge);
                if (n_centre'*ne>0)
                    ne=-ne;
                end
                nx=ne(1);
                ny=ne(2);
                
                
                lx=norm(nodes(elem.nodes(edges.loads(ie,3),1),:)-nodes(elem.nodes(edges.loads(ie,3),2),:));
                ly=norm(nodes(elem.nodes(edges.loads(ie,3),1),:)-nodes(elem.nodes(edges.loads(ie,3),4),:));
                %                 load_Hermite_2D_2
                [p_d1,p_d2,p_d3,p_d4,p_d5,p_d6,p_d7,p_d8,p_d9,p_d10,p_d11,p_d12]=H12_shape_functions(lx,ly);
                
                vx_d1 =-derive_polynom_2D_x_2(p_d1 )/(1j*omega*rho_e);
                vx_d2 =-derive_polynom_2D_x_2(p_d2 )/(1j*omega*rho_e);
                vx_d3 =-derive_polynom_2D_x_2(p_d3 )/(1j*omega*rho_e);
                vx_d4 =-derive_polynom_2D_x_2(p_d4 )/(1j*omega*rho_e);
                vx_d5 =-derive_polynom_2D_x_2(p_d5 )/(1j*omega*rho_e);
                vx_d6 =-derive_polynom_2D_x_2(p_d6 )/(1j*omega*rho_e);
                vx_d7 =-derive_polynom_2D_x_2(p_d7 )/(1j*omega*rho_e);
                vx_d8 =-derive_polynom_2D_x_2(p_d8 )/(1j*omega*rho_e);
                vx_d9 =-derive_polynom_2D_x_2(p_d9 )/(1j*omega*rho_e);
                vx_d10=-derive_polynom_2D_x_2(p_d10)/(1j*omega*rho_e);
                vx_d11=-derive_polynom_2D_x_2(p_d11)/(1j*omega*rho_e);
                vx_d12=-derive_polynom_2D_x_2(p_d12)/(1j*omega*rho_e);
                vy_d1 =-derive_polynom_2D_x_2(p_d1 )/(1j*omega*rho_e);
                vy_d2 =-derive_polynom_2D_x_2(p_d2 )/(1j*omega*rho_e);
                vy_d3 =-derive_polynom_2D_x_2(p_d3 )/(1j*omega*rho_e);
                vy_d4 =-derive_polynom_2D_x_2(p_d4 )/(1j*omega*rho_e);
                vy_d5 =-derive_polynom_2D_x_2(p_d5 )/(1j*omega*rho_e);
                vy_d6 =-derive_polynom_2D_x_2(p_d6 )/(1j*omega*rho_e);
                vy_d7 =-derive_polynom_2D_x_2(p_d7 )/(1j*omega*rho_e);
                vy_d8 =-derive_polynom_2D_x_2(p_d8 )/(1j*omega*rho_e);
                vy_d9 =-derive_polynom_2D_x_2(p_d9 )/(1j*omega*rho_e);
                vy_d10=-derive_polynom_2D_x_2(p_d10)/(1j*omega*rho_e);
                vy_d11=-derive_polynom_2D_x_2(p_d11)/(1j*omega*rho_e);
                vy_d12=-derive_polynom_2D_x_2(p_d12)/(1j*omega*rho_e);
                
                
                
                
                
                
                index_p=dof_A(p_H12(elem.nodes(edges.loads(ie,3),:)));
                indice_test=index_p([1 2 5 4]);
                indice_champs_0=indice_test;
                indice_champs_y=index_p([3 6]);
                
                C=[0 1 0];
                s=1;
                
                
                %                 W=[nx -ny -nx;ny nx -ny;air.Z 0 air.Z];
                %                 Omega=[nx/2 ny/2 1/(2*air.Z);-ny nx 0;-nx/2 -ny/2 1/(2*air.Z)];
                %
                %                 W0plus=W(:,1)
                %                 Wmoins=W(:,3)
                %                 Omega0plus=Omega(1,:);
                %                 Omegamoins=Omega(3,:);
                
                
                
                W_e_in=[nx;ny;air.Z];
                Omega_e_in=[nx/2 ny/2 1/(2*air.Z)];
                W_e_out=[nx;ny;-air.Z];
                Omega_e_out=[nx/2 ny/2 -1/(2*air.Z)];
                
                S_tilde= inv(C*W_e_out)*s;
                R_tilde=-inv(C*W_e_out)*(C*W_e_in);
                
                temp=(W_e_in+W_e_out*R_tilde)*Omega_e_in;
                Boundary_11= (temp(1,:)*nx+temp(2,:)*ny)/(1j*omega);
                temp=W_e_out*S_tilde;
                Boundary_1F= (temp(1,:)*nx+temp(2,:)*ny)/((1j*omega));
                
                
                %                 FF=W_e_out*S_tilde;
                %                 PP=(W_e_in+W_e_out*R_tilde)*Omega_e_in;
                
                
                
                for i_test=1:12
                    eval(['Interp_test=p_d',num2str(i_test),';']);
                    for i_champs=1:12
                        eval(['Interp_champs=Boundary_11(1)*vx_d',num2str(i_champs),'+Boundary_11(2)*vy_d',num2str(i_champs),'+Boundary_11(3)*p_d',num2str(i_champs),';'])
                        A(index_p(i_test),index_p(i_champs))=A(index_p(i_test),index_p(i_champs))-integrate_polynom_2D_edge(multiply_polynom_2D(Interp_test,Interp_champs),a,b);
                    end
                    F(index_p(i_test))=F(index_p(i_test))+integrate_polynom_2D_edge(Interp_test,a,b)*Boundary_1F;
                    
                end
                
                
                %                 for i_test=1:4
                %                     eval(['Psi_test=Psi_',num2str(i_test),'_x;'])
                %                     F(indice_test(i_test))=F(indice_test(i_test))+integrate_polynom(Psi_test,lx_H12)*FF;
                %                     for i_champs=1:4
                %                         eval(['temp_1=PP(3)*Psi_',num2str(i_champs),'_x;']);
                %                         eval(['temp_2=PP(1)*(-1/(j*air.rho*omega))*derive_polynom(Psi_',num2str(i_champs),'_x);']);
                %                         temp=add_polynom(temp_1,temp_2);
                %                         A(indice_test(i_test),indice_champs_0(i_champs))=A(indice_test(i_test),indice_champs_0(i_champs))-integrate_polynom(multiply_polynom(Psi_test,temp),lx_H12);
                %                     end
                %                     for i_champs=1:2
                %                         index_Psi=[1 4];
                %                         %['temp=(1/ly_H12)*PP(2)*(-1/(j*air.rho*omega))*Psi_',num2str(index_Psi(i_champs)),'_x;']
                %                         eval(['temp=(1/ly_H12)*PP(2)*(-1/(j*air.rho*omega))*Psi_',num2str(index_Psi(i_champs)),'_x;'])
                %                         A(indice_test(i_test),indice_champs_y(i_champs))=A(indice_test(i_test),indice_champs_y(i_champs))-integrate_polynom(multiply_polynom(Psi_test,temp),lx_H12);
                %                     end
                %                 end
            else
                aezezaezaezeazezaeazezaeazezaeza
            end
    end
end
