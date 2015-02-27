% msh2DGM.m
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


function [nb_nodes,nb_elements,nb_edges,nodes,elements,element_label,edges,nb_media,num_media,element_num_mat,nb_interfaces,interfaces,nb_MMT,edges_MMT...
    ,nb_loads,loads,nb_dirichlets,dirichlets,nb_periodicity,periodicity,index_label,index_element]=msh2DGM(nom_fichier,tracefigure)

fid=fopen(nom_fichier,'r');
nb_nodes=fscanf(fid,'%i',1);
nb_elements=fscanf(fid,'%i',1);
nb_edges=fscanf(fid,'%i',1);

nodes=zeros(nb_nodes,2);
node_label=zeros(nb_nodes);
elements=zeros(nb_elements,3);
element_label=zeros(nb_elements,1);
edges=zeros(nb_edges,3);


for ii=1:nb_nodes
    nodes(ii,1)=fscanf(fid,'%f',1);
    nodes(ii,2)=fscanf(fid,'%f',1);
    node_label(ii)=fscanf(fid,'%f',1);
end
for ii=1:nb_elements
    elements(ii,1)=fscanf(fid,'%i',1);
    elements(ii,2)=fscanf(fid,'%i',1);
    elements(ii,3)=fscanf(fid,'%i',1);
    element_label(ii,1)=fscanf(fid,'%f',1);
end
for ii=1:nb_edges
    edges(ii,1)=fscanf(fid,'%i',1);
    edges(ii,2)=fscanf(fid,'%i',1);
    edges(ii,3)=fscanf(fid,'%i',1);
end

fclose(fid);

% Creation of 3 segments by element
% segments=[node1 node2 #element1 0 element_label1]
segments=[         elements(:,1) elements(:,2) (1:nb_elements)' zeros(nb_elements,1) element_label];
segments=[segments;elements(:,2) elements(:,3) (1:nb_elements)' zeros(nb_elements,1) element_label];
segments=[segments;elements(:,3) elements(:,1) (1:nb_elements)' zeros(nb_elements,1) element_label];
% Ordering of nodes so that node 1 < node 2
segments(:,1:2)=sort(segments(:,1:2),2);
% Research of double segments
[~, ~, ic] = unique(segments(:,1:2),'rows');
segments=[ic segments(:,1:5)];
% segments=[#segment node1 node2 #element1 0 element_label1]
% Ordering of the vector along the #of segments (first column)
[~,index]=sort(segments(:,1));
segments=segments(index,:);
% Find dupplicate values
temp=find(~diff(segments(:,1)));
% Merging of columns
segments(temp,5)=segments(temp+1,4);
segments(temp,7)=segments(temp+1,6);
% segments=[#segment node1 node2 #element1 #element2 element_label1 element_label2]
segments(temp+1,:)=[];
% Suppression of the first column with temporary index of edges
segments(:,1)=[];
% segments=[node1 node2 #element1 #element2(if any) element_label1 element_label2(if any)]

% Separation between boundary and interfaces;
index_boundary=find(segments(:,4)==0);
index_interface=find(segments(:,4)~=0);

boundaries=[segments(index_boundary,:)];
interfaces=[segments(index_interface,:)];

% check for interfaces that #element1<#element2

i_temp=find(interfaces(:,3)>interfaces(:,4));
temp=interfaces(i_temp,3);
interfaces(i_temp,3)=interfaces(i_temp,4);
interfaces(i_temp,4)=temp;


% interfaces=[node1 node2 #element1 #element2 element_label1 element_label2]
clear segments


% Suppression of temporary values for boundaries
boundaries(:,4:6)=[];
% boundaries=[node1 node2 #element1]

% Ordering of nodes so that column 1 < column 2
edges(:,1:2)=sort(edges(:,1:2),2);

% Merging of boundaries
boundaries=[edges;boundaries];
clear edge
% Research of double segments
[C, ia, ic] = unique(boundaries(:,1:2),'rows');
boundaries=[ic boundaries];
% boundaries=[#boundary node1 node2 #element OR # label]

% Ordering of the vector along the number of boundaries
[temp,index]=sort(boundaries(:,1));
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

% interfaces=[node1 node2 #element1 #element2 0]
nb_interfaces=size(interfaces,1);


temp=unique(element_label(find(floor(element_label/1000)==1)));
nb_media.elas=length(temp);
num_media.elas(1:nb_media.elas)=temp-1000;

temp=unique(element_label(find(floor(element_label/1000)==2)));
nb_media.eqf=length(temp);
num_media.eqf(1:nb_media.eqf)=temp-2000;

temp=unique(element_label(find(floor(element_label/1000)==3)));
nb_media.limp=length(temp);
num_media.limp(1:nb_media.limp)=temp-3000;

temp=unique(element_label(find(floor(element_label/1000)==4)));
nb_media.pem98=length(temp);
num_media.pem98(1:nb_media.pem98)=temp-4000;

temp=unique(element_label(find(floor(element_label/1000)==5)));
nb_media.pem01=length(temp);
num_media.pem01(1:nb_media.pem01)=temp-5000;

temp=unique(element_label(find(element_label==0)));
nb_media.acou=length(temp);


for ie=1:nb_elements
    
    if element_label(ie)==0
        element_num_mat(ie)=0;
    elseif (floor(element_label(ie)/1000)==1)
        element_num_mat(ie)=find(num_media.elas==(element_label(ie)-1000));
        
    elseif (floor(element_label(ie)/1000)==2)
        element_num_mat(ie)=find(num_media.eqf==(element_label(ie)-2000));
    elseif (floor(element_label(ie)/1000)==3)
        element_num_mat(ie)=find(num_media.limp==(element_label(ie)-3000));
        
    elseif (floor(element_label(ie)/1000)==4)
        element_num_mat(ie)=find(num_media.pem98==(element_label(ie)-4000));
    elseif (floor(element_label(ie)/1000)==5)
        element_num_mat(ie)=find(num_media.pem01==(element_label(ie)-5000));
    end
end

% boundaries=[node1 node2 #element #label 0]

temp=find(ismember(abs(boundaries(:,4)),[101]));
edges_MMT=boundaries(temp,:);
boundaries(temp,:)=[];


temp=find(ismember(boundaries(:,4),[1 5 6 9]));
dirichlets=boundaries(temp,:);
boundaries(temp,:)=[];
temp=find(ismember(boundaries(:,4),[98 99]));
periodicity=boundaries(temp,:);
boundaries(temp,:)=[];
temp=find(ismember(boundaries(:,4),[0]));
boundaries(temp,:)=[];


loads=boundaries;


nb_dirichlets=size(dirichlets,1);
nb_loads=size(loads,1);
nb_periodicity=size(periodicity,1);
nb_MMT=size(edges_MMT,1);


index_element=zeros(nb_elements,1);
index_label_temp=zeros(nb_elements,1);
index_element(1)=1;
index_label_temp(1)=element_label(1);
nb_mat=1;
for ielem=2:nb_elements
    for jj=1:(ielem-1)
        bool_test=1;
        if (element_label(ielem)==element_label(jj))
            index_element(ielem)=index_element(jj);
            bool_test=0;
        end
    end
    if (bool_test)
        nb_mat=nb_mat+1;
        index_element(ielem)=nb_mat;
        index_label_temp(nb_mat)=element_label(ielem);
    end
end

index_label=zeros(nb_mat,1);
for ii=1:nb_mat
    index_label(ii)=index_label_temp(ii);
end






%     figure
%     hold on
%
%     for ii=1:nb_loads
%         line([nodes(loads(ii,1),1) nodes(loads(ii,2),1)],[nodes(loads(ii,1),2) nodes(loads(ii,2),2)],'Color','r');
%         line([nodes(loads(ii,2),1) nodes(loads(ii,6),1)],[nodes(loads(ii,2),2) nodes(loads(ii,6),2)],'Color','r');
%         line([nodes(loads(ii,6),1) nodes(loads(ii,1),1)],[nodes(loads(ii,6),2) nodes(loads(ii,1),2)],'Color','r');
%         text(mean(nodes(loads(ii,1:2),1)),mean(nodes(loads(ii,1:2),2)),num2str(loads(ii,3)),'Fontsize',15);
%     end
%         for ii=1:nb_periodicity
%         line([nodes(periodicity(ii,1),1) nodes(periodicity(ii,2),1)],[nodes(periodicity(ii,1),2) nodes(periodicity(ii,2),2)],'Color','r');
%         line([nodes(periodicity(ii,2),1) nodes(periodicity(ii,6),1)],[nodes(periodicity(ii,2),2) nodes(periodicity(ii,6),2)],'Color','r');
%         line([nodes(periodicity(ii,6),1) nodes(periodicity(ii,1),1)],[nodes(periodicity(ii,6),2) nodes(periodicity(ii,1),2)],'Color','r');
%         text(mean(nodes(periodicity(ii,1:2),1)),mean(nodes(periodicity(ii,1:2),2)),num2str(periodicity(ii,3)),'Fontsize',15);
%     end
%
%    axis equal



switch tracefigure
    case 1
        figure
        % Elements
        hold on
        
        for ii=1:nb_elements
            line([nodes(elements(ii,1),1) nodes(elements(ii,2),1)],[nodes(elements(ii,1),2) nodes(elements(ii,2),2)],'Color','r');
            line([nodes(elements(ii,2),1) nodes(elements(ii,3),1)],[nodes(elements(ii,2),2) nodes(elements(ii,3),2)],'Color','r');
            line([nodes(elements(ii,3),1) nodes(elements(ii,1),1)],[nodes(elements(ii,3),2) nodes(elements(ii,1),2)],'Color','r');
            text(mean(nodes(elements(ii,1:3),1)),mean(nodes(elements(ii,1:3),2)),num2str(ii),'Fontsize',15);
        end
        
        plot(nodes(:,1),nodes(:,2),'.','Markersize',15);
        
        
        for ii=1:nb_nodes
            text(nodes(ii,1),nodes(ii,2),num2str(ii),'Fontsize',15);
        end
        axis equal
        
        figure
        hold on
        
        for ii=1:nb_loads
            line([nodes(loads(ii,1),1) nodes(loads(ii,2),1)],[nodes(loads(ii,1),2) nodes(loads(ii,2),2)]);
            text((nodes(loads(ii,1),1)+nodes(loads(ii,2),1))/2,(nodes(loads(ii,1),2)+nodes(loads(ii,2),2))/2,num2str(loads(ii,4)),'Fontsize',15);
        end
        for ii=1:nb_periodicity
            line([nodes(periodicity(ii,1),1) nodes(periodicity(ii,2),1)],[nodes(periodicity(ii,1),2) nodes(periodicity(ii,2),2)]);
            text((nodes(periodicity(ii,1),1)+nodes(periodicity(ii,2),1))/2,(nodes(periodicity(ii,1),2)+nodes(periodicity(ii,2),2))/2,num2str(periodicity(ii,4)),'Fontsize',15);
        end
        for ii=1:nb_MMT
            line([nodes(edges_MMT(ii,1),1) nodes(edges_MMT(ii,2),1)],[nodes(edges_MMT(ii,1),2) nodes(edges_MMT(ii,2),2)]);
            text((nodes(edges_MMT(ii,1),1)+nodes(edges_MMT(ii,2),1))/2,(nodes(edges_MMT(ii,1),2)+nodes(edges_MMT(ii,2),2))/2,num2str(edges_MMT(ii,4)),'Fontsize',15);
        end
        for ii=1:nb_dirichlets
            line([nodes(dirichlets(ii,1),1) nodes(dirichlets(ii,2),1)],[nodes(dirichlets(ii,1),2) nodes(dirichlets(ii,2),2)]);
            text((nodes(dirichlets(ii,1),1)+nodes(dirichlets(ii,2),1))/2,(nodes(dirichlets(ii,1),2)+nodes(dirichlets(ii,2),2))/2,num2str(dirichlets(ii,4)),'Fontsize',15);
        end
        axis equal
        
        
        figure
        hold on
        
        for ii=1:nb_loads
            line([nodes(loads(ii,1),1) nodes(loads(ii,2),1)],[nodes(loads(ii,1),2) nodes(loads(ii,2),2)],'Color','red');
        end
        for ii=1:nb_periodicity
            line([nodes(periodicity(ii,1),1) nodes(periodicity(ii,2),1)],[nodes(periodicity(ii,1),2) nodes(periodicity(ii,2),2)],'Color','k');
            
        end
        for ii=1:nb_MMT
            line([nodes(edges_MMT(ii,1),1) nodes(edges_MMT(ii,2),1)],[nodes(edges_MMT(ii,1),2) nodes(edges_MMT(ii,2),2)]);
            
        end
        for ii=1:nb_interfaces
            line([nodes(interfaces(ii,1),1) nodes(interfaces(ii,2),1)],[nodes(interfaces(ii,1),2) nodes(interfaces(ii,2),2)]);
            
        end
        
        for ii=1:nb_dirichlets
            line([nodes(dirichlets(ii,1),1) nodes(dirichlets(ii,2),1)],[nodes(dirichlets(ii,1),2) nodes(dirichlets(ii,2),2)]);
            
        end
        axis equal
        
        
end

end
