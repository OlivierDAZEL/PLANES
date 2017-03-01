%
% boundary_dispersion_FEM.m
%
% Copyright (C) 2017 Mathieu Gaborit (matael) <mathieu@matael.org>
%
% Licensed under the "THE BEER-WARE LICENSE" (Revision 42):
% Mathieu (matael) Gaborit wrote this file. As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer or coffee in return
%
% Time Convention : exp(+i*omega*t)
%
% 


k_e=omega/air.c;
n_excitation=[cos(data_model.theta_inc); sin(data_model.theta_inc)];
jke=1j*k_e;

e_edge=edges.loads(ie,3);
c_e= mean(nodes(nonzeros(elem.nodes(e_edge,:)),1:2))';
coord_edge(1:2,1)=nodes(edges.loads(ie,1),1:2)';
coord_edge(1:2,2)=nodes(edges.loads(ie,2),1:2)';
a=coord_edge(:,1);
b=coord_edge(:,2);

[nx,ny]=normal_edge_out_element(a,b,c_e);

switch elem.model(e_edge)
    case 1 % TR6
        index_p=dof_A(p_TR(elem.nodes(e_edge,1:6)));
        nb_dof_1=6;
        vcor=nodes(nonzeros(elem.nodes(e_edge,1:6)),1:2);
        base_poly=mean(vcor);
        p_d1=Lagrange_TR6(vcor,1,base_poly);
        p_d2=Lagrange_TR6(vcor,2,base_poly);
        p_d3=Lagrange_TR6(vcor,3,base_poly);
        p_d4=Lagrange_TR6(vcor,4,base_poly);
        p_d5=Lagrange_TR6(vcor,5,base_poly);
        p_d6=Lagrange_TR6(vcor,6,base_poly);
    case 3 % TR3
        not implemented
    case 2
        not implemented
end

factor_rhs = 1j/(omega*air.Z)*(1-[nx, ny]*n_excitation);
factor_lhs = 1j/(omega*air.Z);

for i_test=1:nb_dof_1
    eval(['N_test=p_d',num2str(i_test),';']);

		for i_field=1:nb_dof_1
			eval(['N_field=p_d',num2str(i_field),';']);
			A(index_p(i_field),index_p(i_test))=A(index_p(i_field),index_p(i_test))+factor_lhs*integrate_polynom_2D_edge(N_test,base_poly,N_field,base_poly,a,b,Gauss_points);
		end

    F(index_p(i_test))=F(index_p(i_test))+factor_rhs*integrate_polynom_exp_2D_edge(N_test,base_poly,-jke*n_excitation,[0;0],a,b,Gauss_points);

end
