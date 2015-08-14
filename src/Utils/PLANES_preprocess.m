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

temp=find(ismember(elem.model,[1 3]));
segments=[         elem.nodes(temp,1) elem.nodes(temp,2) temp 0*temp elem.label(temp)];
segments=[segments;elem.nodes(temp,2) elem.nodes(temp,3) temp 0*temp elem.label(temp)];
segments=[segments;elem.nodes(temp,3) elem.nodes(temp,1) temp 0*temp elem.label(temp)];
temp=find(ismember(elem.model,2));
segments=[segments;elem.nodes(temp,1) elem.nodes(temp,2) temp 0*temp elem.label(temp)];
segments=[segments;elem.nodes(temp,2) elem.nodes(temp,3) temp 0*temp elem.label(temp)];
segments=[segments;elem.nodes(temp,3) elem.nodes(temp,4) temp 0*temp elem.label(temp)];
segments=[segments;elem.nodes(temp,4) elem.nodes(temp,1) temp 0*temp elem.label(temp)];


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
temp=find(temp_FEM.*temp_physical);
edges.internal(temp,:)=[];
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

% internal=[node1 node2 #element1 #element2 0]
nb.internal=size(edges.internal,1);


% boundaries=[node1 node2 #element #label 0]

temp=find(ismember(abs(boundaries(:,4)),[101]));
edges.MMT=boundaries(temp,:);
boundaries(temp,:)=[];


temp=find(ismember(boundaries(:,4),[1 5 6 9]));
edges.dirichlets=boundaries(temp,:);
boundaries(temp,:)=[];
temp=find(ismember(boundaries(:,4),[98 99]));
edges.periodicity=boundaries(temp,:);
boundaries(temp,:)=[];
temp=find(ismember(boundaries(:,4),[0]));
boundaries(temp,:)=[];


edges.loads=boundaries;
clear boundaries;

nb.dirichlets=size(edges.dirichlets,1);
nb.loads=size(edges.loads,1);
nb.periodicity=size(edges.periodicity,1);
nb.MMT=size(edges.MMT,1);


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


is_pw=(ismember(edges.loads(:,4),[10 11 12]));
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

is_pw=(ismember(edges.loads(:,4),[13]));
is_pw_R=is_pw;
if sum(is_pw)~=0
    plot_abs=1;
    nb.R=1;
    size_info_vector_R=2;
end


is_pw=(ismember(edges.loads(:,4),[20 21 22]));
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

is_pw=(ismember(edges.loads(:,4),23));
is_pw_T=is_pw;
if sum(is_pw)~=0
    is_pw_T=find(is_pw);
    plot_TL=1;
    nb.T=1;
    size_info_vector_T=2;
end

period=max(nodes(:,1))-min(nodes(:,1));


find_dof_FEM





nb.interfaces=0







if profiles.mesh
   display_mesh 
end



