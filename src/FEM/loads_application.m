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
                [Psi_1_x,Psi_2_x,Psi_3_x,Psi_4_x]=Hermite_shape_functions(lx);                
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
                    F(indice_test(i_test))=F(indice_test(i_test))-integrate_polynom(Psi_test,length_edge)/(1j*omega);
                end
            elseif  ismember(elem.model(edges.loads(ie,3)),[10 11])
                boundary_normal_velocity_fluid
                
            end
            
            
           case {4}  % Unit tangential
            
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
                
                c_1=mean(nodes(nonzeros(edges.loads(ie,4)),1:2))';
                a=nodes(node(1),1:2)';
                b=nodes(node(2),1:2)';
                
                %%%%% vector normal pointing out from e
                
                [nx,ny]=normal_edge_out_element(a,b,c_1);
                tx=-ny;
                ty= nx;
                tx=nx;
                ty=ny;
                

                
                
                F3=TR6_unit(length_edge);
                index_force=dof_A(ux_TR(node));
                
                index_F_elem=find(index_force);
                index_F_global=index_force(index_F_elem);
                F(index_F_global)=F(index_F_global)+F3(index_F_elem)*tx/(1j*omega);
                
                index_force=dof_A(uy_TR(node));
                index_F_elem=find(index_force);
                index_F_global=index_force(index_F_elem);
                F(index_F_global)=F(index_F_global)+F3(index_F_elem)*ty/(1j*omega);
                
      
            end    
            
        case {60}
            if sort(edges.loads(ie,1:2))==sort(elem.nodes(edges.loads(ie,3),1:2)) % Interface node 1 node 2: bottom edge
                stop a revoir
                c_e=air.c;
                k_e=omega/c_e;
                Z_e=air.Z;
                rho_e=air.rho;
                M_e=diag([air.rho,air.rho,1/air.K]);
                
                e_1=edges.loads(ie,3);
                c_1=mean(nodes(nonzeros(elem.nodes(e_1,:)),1:2))';
                
                a=nodes(edges.loads(ie,1),1:2)';
                b=nodes(edges.loads(ie,2),1:2)';
                
                %%%%% vector normal pointing out from e
                
                [nx,ny]=normal_edge_out_element(a,b,c_1);
                
                lx=norm(nodes(elem.nodes(e_1,1),:)-nodes(elem.nodes(e_1,2),:));
                ly=norm(nodes(elem.nodes(e_1,1),:)-nodes(elem.nodes(e_1,4),:));
                %  load_Hermite_2D_2
                [p_d1,p_d2,p_d3,p_d4,p_d5,p_d6,p_d7,p_d8,p_d9,p_d10,p_d11,p_d12]=H12_shape_functions(lx,ly);
                
                vx_d1 =-derive_polynom_2D_x(p_d1 )/(1j*omega*rho_e);
                vx_d2 =-derive_polynom_2D_x(p_d2 )/(1j*omega*rho_e);
                vx_d3 =-derive_polynom_2D_x(p_d3 )/(1j*omega*rho_e);
                vx_d4 =-derive_polynom_2D_x(p_d4 )/(1j*omega*rho_e);
                vx_d5 =-derive_polynom_2D_x(p_d5 )/(1j*omega*rho_e);
                vx_d6 =-derive_polynom_2D_x(p_d6 )/(1j*omega*rho_e);
                vx_d7 =-derive_polynom_2D_x(p_d7 )/(1j*omega*rho_e);
                vx_d8 =-derive_polynom_2D_x(p_d8 )/(1j*omega*rho_e);
                vx_d9 =-derive_polynom_2D_x(p_d9 )/(1j*omega*rho_e);
                vx_d10=-derive_polynom_2D_x(p_d10)/(1j*omega*rho_e);
                vx_d11=-derive_polynom_2D_x(p_d11)/(1j*omega*rho_e);
                vx_d12=-derive_polynom_2D_x(p_d12)/(1j*omega*rho_e);
                vy_d1 =-derive_polynom_2D_x(p_d1 )/(1j*omega*rho_e);
                vy_d2 =-derive_polynom_2D_x(p_d2 )/(1j*omega*rho_e);
                vy_d3 =-derive_polynom_2D_x(p_d3 )/(1j*omega*rho_e);
                vy_d4 =-derive_polynom_2D_x(p_d4 )/(1j*omega*rho_e);
                vy_d5 =-derive_polynom_2D_x(p_d5 )/(1j*omega*rho_e);
                vy_d6 =-derive_polynom_2D_x(p_d6 )/(1j*omega*rho_e);
                vy_d7 =-derive_polynom_2D_x(p_d7 )/(1j*omega*rho_e);
                vy_d8 =-derive_polynom_2D_x(p_d8 )/(1j*omega*rho_e);
                vy_d9 =-derive_polynom_2D_x(p_d9 )/(1j*omega*rho_e);
                vy_d10=-derive_polynom_2D_x(p_d10)/(1j*omega*rho_e);
                vy_d11=-derive_polynom_2D_x(p_d11)/(1j*omega*rho_e);
                vy_d12=-derive_polynom_2D_x(p_d12)/(1j*omega*rho_e);
                
                index_p=dof_A(p_H12(elem.nodes(e_1,:)));
                
                C=[nx ny 0];
                s=-11;
                                
                P_e_in=[nx;ny;air.Z];
                Q_e_in=[nx/2 ny/2 1/(2*air.Z)];
                P_e_out=[nx;ny;-air.Z];
                Q_e_out=[nx/2 ny/2 -1/(2*air.Z)];
                
                S_tilde= inv(C*P_e_out)*s;
                R_tilde=-inv(C*P_e_out)*(C*P_e_in);
                
                temp=(P_e_in+P_e_out*R_tilde)*Q_e_in;
                Boundary_11= (temp(1,:)*nx+temp(2,:)*ny)/(1j*omega);
                temp=P_e_out*S_tilde;
                Boundary_1F= (temp(1,:)*nx+temp(2,:)*ny)/(1j*omega);
                

                for i_test=1:12
                    eval(['Interp_test=p_d',num2str(i_test),';']);
                    for i_champs=1:12
                        eval(['Interp_champs=Boundary_11(1)*vx_d',num2str(i_champs),'+Boundary_11(2)*vy_d',num2str(i_champs),'+Boundary_11(3)*p_d',num2str(i_champs),';'])
                        A(index_p(i_test),index_p(i_champs))=A(index_p(i_test),index_p(i_champs))-integrate_polynom_2D_edge(multiply_polynom_2D(Interp_test,Interp_champs),a,b,Gauss_points);
                    end
                    F(index_p(i_test))=F(index_p(i_test))+integrate_polynom_2D_edge(Interp_test,a,b,Gauss_points)*Boundary_1F;
                    
                end
            else
                aezezaezaezeazezaeazezaeazezaeza
            end
    end
end
