L2=0;

if nb.dof_FEM>0
    L2=L2+X(1:nb.dof_FEM)'*Q_acou*X(1:nb.dof_FEM);
end

if nb.dof_DGM>0
    for ie=1:1:nb.elements
        if ismember(elem.model(ie),[10,11])
            index_DGM=((1:theta_DGM.nb)-1)+dof_start_element(ie);
            c_e=mean(nodes(nonzeros(elem.nodes(ie,:)),1:2))';
            nodes_elem=nodes(nonzeros(elem.nodes(ie,:)),:);
            for i_test=1:theta_DGM.nb
                jktest=+1j*k_e'*[cos(vec_theta(i_test));sin(vec_theta(i_test))];
                for i_field=1:theta_DGM.nb
                    jkfield=-1j*k_e*[cos(vec_theta(i_field));sin(vec_theta(i_field))];
                    L2=L2+(air.Z)^2*conj(X(index_DGM(i_test)))*X(index_DGM(i_field))*integrate_exp_element(nodes_elem,jktest+jkfield,c_e);
                end
            end
            
        end
    end
end
fprintf(file_abs_id,'%d\t%1.15e\t%1.15e\n',nb.dof_FEM+nb.dof_DGM,frequency.vec(i_f),L2);
