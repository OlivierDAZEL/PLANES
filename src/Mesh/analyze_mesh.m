is_pw=(ismember(loads(:,4),[10 11 12]));
is_pw_R=is_pw;
if sum(is_pw)~=0
    plot_abs=1;
    nb_R=1;
    size_info_vector_R=1;
else
    plot_abs=0;
    nb_R=0;
    size_info_vector_R=1;
end

is_pw=(ismember(loads(:,4),[13]));
is_pw_R=is_pw;
if sum(is_pw)~=0
    plot_abs=1;
    nb_R=1;
    size_info_vector_R=2;
end


is_pw=(ismember(loads(:,4),[20 21 22]));
is_pw_T=is_pw;
if sum(is_pw)~=0
    is_pw_T=find(is_pw);
    plot_TL=1;
    nb_T=1;
    size_info_vector_T=1;
else
    plot_TL=0;
    nb_T=0;
    size_info_vector_T=1;
end

is_pw=(ismember(loads(:,4),23));
is_pw_T=is_pw;
if sum(is_pw)~=0
    is_pw_T=find(is_pw);
    plot_TL=1;
    nb_T=1;
    size_info_vector_T=2;
end



period=max(nodes(:,1))-min(nodes(:,1));

check_dirichlet



if exist('incident')
    for ii=1:length(incident)
        if floor(incident(ii).typ/1000)==1
           eval(['Mat_elas_' num2str(incident(1).typ-1000)])
           incident(ii).rho=rho_solide;
           incident(ii).lambda=lambda_solide;
           incident(ii).mu=mu_solide;
        end
            
    end
end

if exist('transmission')
    for ii=1:length(transmission)
        if floor(transmission(ii).typ/1000)==1
           eval(['Mat_elas_' num2str(transmission(1).typ-1000)])
           transmission(ii).rho=rho_solide;
           transmission(ii).lambda=lambda_solide;
           transmission(ii).mu=mu_solide;
        end
            
    end
end

%node_MMT=[];
if nb_MMT~=0
    
    index_MMT_moins=find(edges_MMT(:,4)<0);
    
    edges_MMT_moins=edges_MMT(index_MMT_moins,:);
    for ii=1:length(index_MMT_moins);
        
        % Association of boundaries on plus and on minus by the middle node
        node_moins=edges_MMT(index_MMT_moins(ii),6);
        [~,node_plus]=min(abs((nodes(:,1)-nodes(node_moins,1)-delta_x_MMT)+1i*(nodes(:,2)-nodes(node_moins,2)-delta_y_MMT)));
        
        index_MMT_plus=find(edges_MMT(:,6)==node_plus);
        
        edges_MMT_plus(ii,:)=edges_MMT(index_MMT_plus,:);
        %node_MMT=[node_MMT; node2_moins node2_plus];
        
        
        node_moins=edges_MMT(index_MMT_moins(ii),1);
        [~,node_plus]=min(abs((nodes(:,1)-nodes(node_moins,1)-delta_x_MMT)+1i*(nodes(:,2)-nodes(node_moins,2)-delta_y_MMT)));

        node1_plus=find(edges_MMT_plus(ii,1:2)==node_plus);
        
        if (~isempty(node1_plus))
            if node1_plus==2
               edges_MMT_plus(ii,1:2)=edges_MMT_plus(ii,[2 1]); 
            end
            
        else
             stop
        end
    end
    
    temp=element_label(edges_MMT_moins(:,3));
    temp=temp-temp(1)*ones(size(temp));
    if norm(temp)==0
       element_MMT_moins=element_label(edges_MMT_moins(1,3));
    else
        stop
    end
    
    temp=element_label(edges_MMT_plus(:,3));
    temp=temp-temp(1)*ones(size(temp));
    if norm(temp)==0
       element_MMT_plus=element_label(edges_MMT_plus(1,3));
    else
        stop
    end
    


    

end

