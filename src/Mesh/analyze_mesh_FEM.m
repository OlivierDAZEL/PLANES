% analyze_mesh_FEM.m
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

if nb_periodicity~=0
    
    edge_left= find(periodicity(:,4)==98);
    edge_right=find(periodicity(:,4)==99);
    
    node_left=unique([periodicity(edge_left,1);periodicity(edge_left,2);periodicity(edge_left,6)]);
    [temp,i_left]=sort(nodes(node_left,2));
    node_left=node_left(i_left);
    
    node_right=unique([periodicity(edge_right,1);periodicity(edge_right,2);periodicity(edge_right,6)]);
    [temp,i_right]=sort(nodes(node_right,2));
    node_right=node_right(i_right);
    
    
    dof_left= dof_A([3*(node_left -1)+1;3*(node_left -1)+2;3*(node_left -1)+3]);
    dof_right=dof_A([3*(node_right-1)+1;3*(node_right-1)+2;3*(node_right-1)+3]);
    
    dof_left= dof_left (find(dof_left));
    dof_right=dof_right(find(dof_right));
    
end


