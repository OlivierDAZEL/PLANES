



% Computation of the reference solution
B=-1/(1j*omega*sin(-k_air*model_data.ly));
%p_air=-B_anal*(air.K)*k_air*cos(k_air*x_air);
X_ref=zeros(nb.dof_FEM+nb.dof_DGM,1);
nb_theta_ref=2;
vec_theta_ref=[pi/2 -pi/2;];



for ie=1:nb.elements
    if ismember(elem.model(ie),[1 2 3])
        y=nodes(nonzeros(elem.nodes(ie,:)),2);
        p_ref=-B*(air.K)*k_air*cos(k_air*(y-model_data.ly));
        if elem.model(ie)==2
            ly=norm(nodes(elem.nodes(ie,1),:)-nodes(elem.nodes(ie,4),:));
            index_p=dof_A(p_H12(elem.nodes(ie,1:4)));
            X_ref(index_p(1:3:end))=p_ref;
            dp_ref_dy=ly*B*(air.K)*k_air^2*sin(k_air*(y-model_data.ly));
            X_ref(index_p(3:3:end))=dp_ref_dy;
        else
            index_p=dof_A(p_TR(nonzeros(elem.nodes(ie,:))));
            X_ref(index_p)=p_ref;
        end
    end
    if ismember(elem.model(ie),[10 11])
        center_element=mean(nodes(nonzeros(elem.nodes(ie,:)),1:2))';
        X_ref(dof_start_element(ie)             )  =(-B*(air.K)*k_air*exp( 1j*k_air*model_data.ly)/2)*exp(-1j*k_air*center_element(2))/air.Z;
        X_ref(dof_start_element(ie)+1)=(-B*(air.K)*k_air*exp(-1j*k_air*model_data.ly)/2)*exp( 1j*k_air*center_element(2))/air.Z;
    end
end
L2_ref=((model_data.lx*(air.K*k_air*abs(B))^2)/2)*(model_data.ly+sin(2*k_air*model_data.ly)/(2*k_air));

sol_ref(dof_back)=X_ref(1:nb.dof_FEM);


L2=0;
L2_error=0;

if nb.dof_FEM>0
    L2=L2+X(1:nb.dof_FEM)'*Q_acou*X(1:nb.dof_FEM);
    L2_error=L2_error+(X(1:nb.dof_FEM)-X_ref(1:nb.dof_FEM))'*Q_acou*(X(1:nb.dof_FEM)-X_ref(1:nb.dof_FEM));
end


for ie=1:1:nb.elements
    if ismember(elem.model(ie),[10,11])
        index_DGM=((1:theta_DGM.nb)-1)+dof_start_element(ie);
        c_e=mean(nodes(nonzeros(elem.nodes(ie,:)),1:2))';
        nodes_elem=nodes(nonzeros(elem.nodes(ie,:)),:);
        for i_test=1:theta_DGM.nb
            jktest=+1j*k_e'*[cos(vec_theta(i_test));sin(vec_theta(i_test))];
            for i_field=1:theta_DGM.nb
                jkfield=-1j*k_e*[cos(vec_theta(i_field));sin(vec_theta(i_field))];
                L2_error=L2_error+(air.Z)^2*conj(X(index_DGM(i_test)))*X(index_DGM(i_field))*integrate_exp_element(nodes_elem,jktest+jkfield,c_e);
                L2=L2+(air.Z)^2*conj(X(index_DGM(i_test)))*X(index_DGM(i_field))*integrate_exp_element(nodes_elem,jktest+jkfield,c_e);
            end
            for i_field=1:nb_theta_ref
                jkfield=-1j*k_e*[cos(vec_theta_ref(i_field));sin(vec_theta_ref(i_field))];
                L2_error=L2_error-(air.Z)^2*conj(X(index_DGM(i_test)))*X_ref(index_DGM(i_field))*integrate_exp_element(nodes_elem,jktest+jkfield,c_e);
            end
        end
        for i_test=1:nb_theta_ref
            jktest=+1j*k_e'*[cos(vec_theta_ref(i_test));sin(vec_theta_ref(i_test))];
            for i_field=1:theta_DGM.nb
                jkfield=-1j*k_e*[cos(vec_theta(i_field));sin(vec_theta(i_field))];
                L2_error=L2_error-(air.Z)^2*conj(X_ref(index_DGM(i_test)))*X(index_DGM(i_field))*integrate_exp_element(nodes_elem,jktest+jkfield,c_e);
            end
            for i_field=1:nb_theta_ref
                jkfield=-1j*k_e*[cos(vec_theta_ref(i_field));sin(vec_theta_ref(i_field))];
                L2_error=L2_error+(air.Z)^2*conj(X_ref(index_DGM(i_test)))*X_ref(index_DGM(i_field))*integrate_exp_element(nodes_elem,jktest+jkfield,c_e);
            end
        end
        
    end
end