% PLANES_preprocess.m
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

%% First step: identify boudaries and internal edges
% Creation of 3 segments by element
% segments=[node1 node2 #element1 0 elem.label1]

init_vec_frequencies



if exist('nodes')
    segments=zeros(1,5); % Line 1 to be removed at the end
    
    temp=find(ismember(elem.model,[1 3 10]));
    temp=reshape(temp,length(temp),1);
    if length(temp)>0
        segments=[segments;elem.nodes(temp,1) elem.nodes(temp,2) temp 0*temp elem.label(temp)];
        segments=[segments;elem.nodes(temp,2) elem.nodes(temp,3) temp 0*temp elem.label(temp)];
        segments=[segments;elem.nodes(temp,3) elem.nodes(temp,1) temp 0*temp elem.label(temp)];
    end
    temp=find(ismember(elem.model,[2 11]));
    temp=reshape(temp,length(temp),1);
    
    if length(temp)>0
        segments=[segments;elem.nodes(temp,1) elem.nodes(temp,2) temp 0*temp elem.label(temp)];
        segments=[segments;elem.nodes(temp,2) elem.nodes(temp,3) temp 0*temp elem.label(temp)];
        segments=[segments;elem.nodes(temp,3) elem.nodes(temp,4) temp 0*temp elem.label(temp)];
        segments=[segments;elem.nodes(temp,4) elem.nodes(temp,1) temp 0*temp elem.label(temp)];
    end
    segments(1,:)=[]; % Line 1 is removed at the end
    
    % Ordering of nodes so that node 1 < node 2
    segments(:,1:2)=sort(segments(:,1:2),2);
    % Research of double segments
    [~, ~, temp] = unique(segments(:,1:2),'rows');
    segments=[temp segments(:,1:5)];
    % segments=[#segment node1 node2 #element1 0 elem.label1]
    % Ordering of the vector along the #of segments (first column)
    [~,temp]=sort(segments(:,1));
    segments=segments(temp,:);
    % Find dupplicate values
    temp=find(~diff(segments(:,1)));
    % Merging of columns
    segments(temp,5)=segments(temp+1,4);
    segments(temp,7)=segments(temp+1,6);
    % segments=[#segment node1 node2 #element1 #element2 elem.label1 elem.label2]
    segments(temp+1,:)=[];
    % Suppression of the first column with temporary index of edges
    segments(:,1)=[];
    % segments=[node1 node2 #element1 #element2(if any) elem.label1 elem.label2(if any)]
    
    % Separation between boundary and internal;
    edges.internal=[segments(find(segments(:,4)~=0),:)];
    % check for internal that #element1<#element2
    temp=find(edges.internal(:,3)>edges.internal(:,4));
    v_temp=edges.internal(temp,3);
    edges.internal(temp,3)=edges.internal(temp,4);
    edges.internal(temp,4)=v_temp;
    
    % Check if the elements are both DGM
    temp=find(ismember(elem.model(edges.internal(:,3)),[10 11]).*ismember(elem.model(edges.internal(:,4)),[10 11]));
    
    edges.internal_DGM=edges.internal(temp,1:4);
    nb.internal_DGM=length(temp);
    
    edges.internal(temp,:)=[];
    
    
    boundaries=    [segments(find(segments(:,4)==0),:)];
    
    
    clear segments
    
    % internal=[node1 node2 #element1 #element2 elem.label1 elem.label2]
    
    
    %% Second step: for internal edges : identify for internal edges the boundary with natural coupling
    %%% They correspond to internal edges with FEM in both elements on naturally coupling physical medium
    
    temp_FEM=(ismember(elem.model(edges.internal(:,3)),1).*ismember(elem.model(edges.internal(:,4)),1));
    temp_FEM=temp_FEM+(ismember(elem.model(edges.internal(:,3)),2).*ismember(elem.model(edges.internal(:,4)),2));
    temp_FEM=temp_FEM+(ismember(elem.model(edges.internal(:,3)),3).*ismember(elem.model(edges.internal(:,4)),3));
    %The media on both faces of the internal edge are modelled both by FEM
    % Check of they are of same physical nature
    temp_physical=(ismember(floor(elem.label(edges.internal(:,3))/1000),[0 2 3])).*(ismember(floor(elem.label(edges.internal(:,4))/1000),[0 2 3]));
    temp_physical=temp_physical+(ismember(floor(elem.label(edges.internal(:,3))/1000),[1])).*(ismember(floor(elem.label(edges.internal(:,4))/1000),[1]));
    temp_physical=temp_physical+(ismember(floor(elem.label(edges.internal(:,3))/1000),[4])).*(ismember(floor(elem.label(edges.internal(:,4))/1000),[4]));
    temp_physical=temp_physical+(ismember(floor(elem.label(edges.internal(:,3))/1000),[5])).*(ismember(floor(elem.label(edges.internal(:,4))/1000),[5]));
    temp_physical=temp_physical+(ismember(floor(elem.label(edges.internal(:,3))/1000),[1])).*(ismember(floor(elem.label(edges.internal(:,4))/1000),[5]));
    temp_physical=temp_physical+(ismember(floor(elem.label(edges.internal(:,3))/1000),[5])).*(ismember(floor(elem.label(edges.internal(:,4))/1000),[1]));
    temp=find(temp_FEM.*temp_physical);
    edges.internal(temp,:)=[];
    
    % internal=[node1 node2 #element1 #element2 0]
    nb.internal=size(edges.internal,1);
    
    clear temp_physical
    
    
    
    % Suppression of temporary values for boundaries
    boundaries(:,4:6)=[];
    % boundaries=[node1 node2 #element1]
    
    
    % Ordering of nodes so that column 1 < column 2
    edge_msh(:,1:2)=sort(edge_msh(:,1:2),2);
    
    % Merging of boundaries
    boundaries=[edge_msh;boundaries];
    
    
    clear edge_msh
    % Research of double segments
    [~,~,temp] = unique(boundaries(:,1:2),'rows');
    boundaries=[temp boundaries];
    
    
    % boundaries=[#boundary node1 node2 #element OR # label]
    
    % Ordering of the vector along the number of boundaries
    [~,index]=sort(boundaries(:,1));
    boundaries=boundaries(index,:);
    
    % Find dupplicate values
    temp=find(~diff(boundaries(:,1)));
    % Merging of columns
    boundaries(temp,5)=boundaries(temp+1,4);
    % boundaries=[#boundary node1 node2 #label #element]
    boundaries(temp+1,:)=[];
    % Suppression of the first column with temporary index of boundaries
    boundaries(:,1)=[];
    
    % boundaries=[node1 node2 #label #element]
    boundaries(:,[4 3])=boundaries(:,[3 4]);
    % boundaries=[node1 node2 #element #label]
    
    
    
    temp=find(ismember(abs(boundaries(:,4)),[0]));
    
    boundaries(temp,:)=[];
    
    
    % Suppression of interfaces with FEM and natural coupling
    
    
    % Find Dirichlet boundaries : label 1 not FEM or 5 6 9
    
    temp=      (ismember(boundaries(:,4),[1]));
    temp=temp.*ismember(elem.model(boundaries(:,3)),[1 2 3]);
    temp=find(temp);
    boundaries(temp,:)=[];
    
    
    
    temp=      (ismember(boundaries(:,4),[1 5 6 9]));
    edges.dirichlets=boundaries(temp,:);
    boundaries(temp,:)=[];
    
    
    
    temp=find(ismember(boundaries(:,4),[25]));
    edges.incompatible=boundaries(temp,:);
    
    
    boundaries(temp,:)=[];
    
    
    identify_incompatible_mesh
    
    edges.flux=[edges.flux;edges.internal_DGM];
    nb.flux=size(edges.flux,1);
    nb.internal_DGM=0;
    edges=rmfield(edges,'internal_DGM');
    
    
    temp=find(ismember(boundaries(:,4),[98 99]));
    edges.periodicity=boundaries(temp,:);
    boundaries(temp,:)=[];
    
    temp=find(ismember(boundaries(:,4),[10 11 12 13 20 21 22 23]));
    edges.DtN=boundaries(temp,:);
    boundaries(temp,:)=[];
    
    
    temp=find(ismember(floor(boundaries(:,4)/100),[4]));
    edges.ZOD=boundaries(temp,:);
    boundaries(temp,:)=[];
    
    
    temp=find(ismember(boundaries(:,4),[0]));
    boundaries(temp,:)=[];
    
    
    
    edges.loads=boundaries;
    clear boundaries;
    
    nb.dirichlets=size(edges.dirichlets,1);
    nb.loads=size(edges.loads,1);
    nb.periodicity=size(edges.periodicity,1);
    nb.ZOD=size(edges.ZOD,1);
    nb.DtN=size(edges.DtN,1);
    
    
    if (sum(elem.model==1)~=0)
        [nb,nodes,elem,edges]=TR32TR6(nb,nodes,elem,edges);
    end
    
    if (sum(elem.model==2)~=0)
        [elem,H_elem_H12,Q_elem_H12]=create_elementary_H12(nb,nodes,elem);
    end
    
    
    
    temp=unique(elem.label(find(floor(elem.label/1000)==1)));
    nb.media.elas=length(temp);
    num_media.elas(1:nb.media.elas)=temp-1000;
    
    temp=unique(elem.label(find(floor(elem.label/1000)==2)));
    nb.media.eqf=length(temp);
    num_media.eqf(1:nb.media.eqf)=temp-2000;
    
    temp=unique(elem.label(find(floor(elem.label/1000)==3)));
    nb.media.limp=length(temp);
    num_media.limp(1:nb.media.limp)=temp-3000;
    
    temp=unique(elem.label(find(floor(elem.label/1000)==4)));
    nb.media.pem98=length(temp);
    num_media.pem98(1:nb.media.pem98)=temp-4000;
    
    temp=unique(elem.label(find(floor(elem.label/1000)==5)));
    nb.media.pem01=length(temp);
    num_media.pem01(1:nb.media.pem01)=temp-5000;
    
    temp=unique(elem.label(find(floor(elem.label/1000)==0)));
    nb.media.acou=length(temp);
    
    temp=unique(elem.label(find(floor(elem.label/1000)==8)));
    nb.media.PML=length(temp);
    
    for ie=1:nb.elements
        
        if (floor(elem.label(ie)/1000)==0)
            elem.num_mat(ie)=0;
        elseif (floor(elem.label(ie)/1000)==1)
            elem.num_mat(ie)=find(num_media.elas==(elem.label(ie)-1000));
        elseif (floor(elem.label(ie)/1000)==2)
            elem.num_mat(ie)=find(num_media.eqf==(elem.label(ie)-2000));
        elseif (floor(elem.label(ie)/1000)==3)
            elem.num_mat(ie)=find(num_media.limp==(elem.label(ie)-3000));
            
        elseif (floor(elem.label(ie)/1000)==4)
            elem.num_mat(ie)=find(num_media.pem98==(elem.label(ie)-4000));
        elseif (floor(elem.label(ie)/1000)==5)
            elem.num_mat(ie)=find(num_media.pem01==(elem.label(ie)-5000));
        end
    end
    
    
    is_pw=(ismember(edges.DtN(:,4),[10 11 12]));
    is_pw_R=is_pw;
    if sum(is_pw)~=0
        plot_abs=1;
        nb.R=1;
        size_info_vector_R=1;
    else
        plot_abs=0;
        nb.R=0;
        size_info_vector_R=1;
    end
    
    is_pw=(ismember(edges.DtN(:,4),[13]));
    is_pw_R=is_pw;
    if sum(is_pw)~=0
        plot_abs=1;
        nb.R=1;
        size_info_vector_R=2;
    end
    
    
    is_pw=(ismember(edges.DtN(:,4),[20 21 22]));
    is_pw_T=is_pw;
    if sum(is_pw)~=0
        is_pw_T=find(is_pw);
        plot_TL=1;
        nb.T=1;
        size_info_vector_T=1;
    else
        plot_TL=0;
        nb.T=0;
        size_info_vector_T=1;
    end
    
    is_pw=(ismember(edges.DtN(:,4),23));
    is_pw_T=is_pw;
    if sum(is_pw)~=0
        is_pw_T=find(is_pw);
        plot_TL=1;
        nb.T=1;
        size_info_vector_T=2;
    end
    
    
    if nb.ZOD~=0
        
        index_ZOD_moins=find(mod(edges.ZOD(:,4),2)==1);
        
        edges.ZOD_moins=edges.ZOD(index_ZOD_moins,:);
        
        
        
        for ii=1:length(index_ZOD_moins);
            number_ZOD=1+(edges.ZOD(index_ZOD_moins(ii),4)-401)/2;
            % Association of boundaries on plus and on minus by the middle node
            node_moins=edges.ZOD(index_ZOD_moins(ii),6);
            [~,node_plus]=min(abs((nodes(:,1)-nodes(node_moins,1)-data_model.multilayer_ZOD(number_ZOD).delta_x)+1i*(nodes(:,2)-nodes(node_moins,2)-data_model.multilayer_ZOD(number_ZOD).delta_y)));
            
            index_ZOD_plus=find(edges.ZOD(:,6)==node_plus);
            
            edges.ZOD_plus(ii,:)=edges.ZOD(index_ZOD_plus,:);
            
            
            node_moins=edges.ZOD(index_ZOD_moins(ii),1);
            [~,node_plus]=min(abs((nodes(:,1)-nodes(node_moins,1)-data_model.multilayer_ZOD(number_ZOD).delta_x)+1i*(nodes(:,2)-nodes(node_moins,2)-data_model.multilayer_ZOD(number_ZOD).delta_y)));
            
            node1_plus=find(edges.ZOD_plus(ii,1:2)==node_plus);
            
            if (~isempty(node1_plus))
                if node1_plus==2
                    edges.ZOD_plus(ii,1:2)=edges.ZOD_plus(ii,[2 1]);
                end
                
            else
                stop
            end
        end
        
        %         temp=elem.label(edges.ZOD_moins(:,3));
        %         temp=temp-temp(1)*ones(size(temp));
        %         if norm(temp)==0
        %             elem.ZOD_moins=elem.label(edges.ZOD_moins(1,3));
        %         else
        %             stop
        %         end
        %
        %         temp=elem.label(edges.ZOD_plus(:,3));
        %         temp=temp-temp(1)*ones(size(temp));
        %         if norm(temp)==0
        %             elem.ZOD_plus=elem.label(edges.ZOD_plus(1,3));
        %         else
        %             stop
        %         end
    end
    
    period=max(nodes(:,1))-min(nodes(:,1));
    
    
    find_dof_FEM
    find_dof_DGM
    
    I_inc=zeros(frequency.nb,1);
    W_vis=zeros(frequency.nb,1);
    W_struct=zeros(frequency.nb,1);
    W_therm=zeros(frequency.nb,1);
    W_elas=zeros(frequency.nb,1);
    abs_vis=zeros(frequency.nb,1);
    abs_struct=zeros(frequency.nb,1);
    abs_therm=zeros(frequency.nb,1);
    abs_elas=zeros(frequency.nb,1);
    TL_EF=zeros(frequency.nb,1);
    abs_EF=zeros(frequency.nb,1);
    
    if isfield(data_model,'theta_DGM')
        vec_theta=data_model.tilt+linspace(0,2*pi,data_model.theta_DGM.nb+1);
        vec_theta(end)=[];
        vec_theta=vec_theta+pi/2;
    end
    
    
    if data_model.profiles.mesh
        display_mesh
    end
    
    if (nb.R~=0)
        file_abs_id=fopen(name.file_abs,'w');
    end
    if (nb.T~=0)
        file_TL_id=fopen(name.file_TL,'w');
    end
else
    nb.dof_FEM=0;
    nb.dof_DGM=0;
end

if exist('multilayer')
    nb_multilayers=size(multilayer,2);
    % Addition of a new layer for the incident medium
    l0=multilayer(1,:);
    for i_m=1:nb_multilayers
        l0(1,i_m).mat=0;
        l0(1,i_m).d=0;
        l0(1,i_m).nb=l0(1,i_m).nb+1;
    end
    multilayer=[l0;multilayer];
end



