% createmshH12DGM.m
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


function [nb,nodes,elements,element_label,element_model,edges,num_media,element_num_mat,interfaces,edges_MMT,loads,dirichlets,periodicity]=createmshH16(lx,ly,nx,ny,tracefigure)

nb.nodes=(nx+1)*(ny+1);
nb.elements=nx*ny;
nb.edges=2*(nx+ny);

nodes=zeros(nb.nodes,2);
node_label=zeros(nb.nodes);
elements=zeros(nb.elements,3);
element_label=zeros(nb.elements,1);
element_model=zeros(nb.elements,1);
edges=zeros(nb.edges,3);

temp=0;
Delta_x=lx/nx;
Delta_y=ly/ny;
for iy=1:ny+1
    ynode=(iy-1)*Delta_y;
    
    for ix=1:nx+1
        xnode=(ix-1)*Delta_x;
        temp=temp+1;
        nodes(temp,1)=xnode;
        nodes(temp,2)=ynode;
        node_label(temp)=0;
    end
end

temp=0;
for ix=1:nx
    for iy=1:ny
        temp=temp+1;
        node_base=ix+(iy-1)*(nx+1);
        elements(temp,1)=node_base;
        elements(temp,2)=node_base+1;
        elements(temp,3)=node_base+(nx+2);
        elements(temp,4)=node_base+(nx+1);
        element_label(temp)=0;
    end
end




temp=0;
for ix=1:nx
    temp=temp+1;
    node_base=ix;
    edges(temp,1)=node_base;
    edges(temp,2)=node_base+1;
    edges(temp,3)=3;
end

for ix=1:nx
    temp=temp+1;
    node_base=ix+ny*(nx+1);
    edges(temp,1)=node_base;
    edges(temp,2)=node_base+1;
    edges(temp,3)=1;
end

for iy=1:ny
    temp=temp+1;
    node_base=nx+1+(iy-1)*(nx+1);
    edges(temp,1)=node_base;
    edges(temp,2)=node_base+nx+1;
    edges(temp,3)=1;
end

for iy=1:ny
    temp=temp+1;
    node_base=1+(iy-1)*(nx+1);
    edges(temp,1)=node_base;
    edges(temp,2)=node_base+nx+1;
    edges(temp,3)=1;
end



% Creation of 4 segments by element
% segments=[node1 node2 #element1 0 element_label1]
segments=[         elements(:,1) elements(:,2) (1:nb.elements)' zeros(nb.elements,1) element_label];
segments=[segments;elements(:,2) elements(:,3) (1:nb.elements)' zeros(nb.elements,1) element_label];
segments=[segments;elements(:,3) elements(:,4) (1:nb.elements)' zeros(nb.elements,1) element_label];
segments=[segments;elements(:,4) elements(:,1) (1:nb.elements)' zeros(nb.elements,1) element_label];
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




% interfaces=[node1 node2 #element1 #element2 element_label1 element_label2]
clear segments
% Research of interfaces between two similar elements
temp=find(floor(interfaces(:,5)/1000)==floor(interfaces(:,6)/1000));
%Suppressions of these interfaces
interfaces(temp,:)=[];
% Research of interfaces between PML and air
temp=find((interfaces(:,5)==0)&(floor(interfaces(:,6)/1000)==8)|(interfaces(:,6)==0)&(floor(interfaces(:,5)/1000)==8));
%Suppressions of these interfaces
interfaces(temp,:)=[];
% Research of interfaces between air and Biot 1998
temp=find((interfaces(:,5)==0)&(floor(interfaces(:,6)/1000)==4)|(interfaces(:,6)==0)&(floor(interfaces(:,5)/1000)==4));
%Suppressions of these interfaces
interfaces(temp,:)=[]
nb.interfaces=size(interfaces,1);


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




temp=unique(element_label(find(floor(element_label/1000)==1)));
nb.media.elas=length(temp);
num_media.elas(1:nb.media.elas)=temp-1000;

temp=unique(element_label(find(floor(element_label/1000)==2)));
nb.media.eqf=length(temp);
num_media.eqf(1:nb.media.eqf)=temp-2000;

temp=unique(element_label(find(floor(element_label/1000)==3)));
nb.media.limp=length(temp);
num_media.limp(1:nb.media.limp)=temp-3000;

temp=unique(element_label(find(floor(element_label/1000)==4)));
nb.media.pem98=length(temp);
num_media.pem98(1:nb.media.pem98)=temp-4000;

temp=unique(element_label(find(floor(element_label/1000)==5)));
nb.media.pem01=length(temp);
num_media.pem01(1:nb.media.pem01)=temp-5000;

temp=unique(element_label(find(floor(element_label/1000)==0)));
nb.media.acou=length(temp);

temp=unique(element_label(find(floor(element_label/1000)==8)));
nb.media.PML=length(temp);

for ie=1:nb.elements
    if (floor(element_label(ie)/1000)==0)
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

% boundaries=[node1 node2 #element #label 0 node_middle]



temp=find(ismember(abs(boundaries(:,4)),[101]));
edges_MMT=boundaries(temp,:);
boundaries(temp,:)=[];


temp=find(ismember(boundaries(:,4),[1 5 6]));
dirichlets=boundaries(temp,:);
boundaries(temp,:)=[];
temp=find(ismember(boundaries(:,4),[98 99]));
periodicity=boundaries(temp,:);
boundaries(temp,:)=[];
temp=find(ismember(boundaries(:,4),[0]));
boundaries(temp,:)=[];


loads=boundaries;


nb.dirichlets=size(dirichlets,1);
nb.loads=size(loads,1);
nb.periodicity=size(periodicity,1);
nb.MMT=size(edges_MMT,1);

%     figure
%     hold on
%
%     for ii=1:nb.loads
%         line([nodes(loads(ii,1),1) nodes(loads(ii,2),1)],[nodes(loads(ii,1),2) nodes(loads(ii,2),2)],'Color','r');
%         line([nodes(loads(ii,2),1) nodes(loads(ii,6),1)],[nodes(loads(ii,2),2) nodes(loads(ii,6),2)],'Color','r');
%         line([nodes(loads(ii,6),1) nodes(loads(ii,1),1)],[nodes(loads(ii,6),2) nodes(loads(ii,1),2)],'Color','r');
%         text(mean(nodes(loads(ii,1:2),1)),mean(nodes(loads(ii,1:2),2)),num2str(loads(ii,3)),'Fontsize',15);
%     end
%         for ii=1:nb.periodicity
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
        hold on
        
        for ii=1:nb.elements
            line([nodes(elements(ii,1),1) nodes(elements(ii,2),1)],[nodes(elements(ii,1),2) nodes(elements(ii,2),2)],'Color','r');
            line([nodes(elements(ii,2),1) nodes(elements(ii,3),1)],[nodes(elements(ii,2),2) nodes(elements(ii,3),2)],'Color','r');
            line([nodes(elements(ii,3),1) nodes(elements(ii,4),1)],[nodes(elements(ii,3),2) nodes(elements(ii,4),2)],'Color','r');
            line([nodes(elements(ii,4),1) nodes(elements(ii,1),1)],[nodes(elements(ii,4),2) nodes(elements(ii,1),2)],'Color','r');
            text(mean(nodes(elements(ii,:),1)),mean(nodes(elements(ii,:),2)),num2str(element_label(ii)),'Fontsize',15);
        end
        
        plot(nodes(:,1),nodes(:,2),'.','Markersize',15);
        for ii=1:nb.nodes
            text(nodes(ii,1),nodes(ii,2),num2str(ii),'Fontsize',15);
        end
        
        
        %     for ii=1:nb.nodes
        %         text(nodes(ii,1),nodes(ii,2),num2str(ii),'Fontsize',15);
        %     end
        axis equal
        
        figure
        hold on
        
        for ii=1:nb.loads
            line([nodes(loads(ii,1),1) nodes(loads(ii,2),1)],[nodes(loads(ii,1),2) nodes(loads(ii,2),2)]);
            text((nodes(loads(ii,1),1)+nodes(loads(ii,2),1))/2,(nodes(loads(ii,1),2)+nodes(loads(ii,2),2))/2,num2str(loads(ii,4)),'Fontsize',15);
        end
        for ii=1:nb.periodicity
            line([nodes(periodicity(ii,1),1) nodes(periodicity(ii,2),1)],[nodes(periodicity(ii,1),2) nodes(periodicity(ii,2),2)]);
            text((nodes(periodicity(ii,1),1)+nodes(periodicity(ii,2),1))/2,(nodes(periodicity(ii,1),2)+nodes(periodicity(ii,2),2))/2,num2str(periodicity(ii,4)),'Fontsize',15);
            
        end
        for ii=1:nb.MMT
            line([nodes(edges_MMT(ii,1),1) nodes(edges_MMT(ii,2),1)],[nodes(edges_MMT(ii,1),2) nodes(edges_MMT(ii,2),2)]);
            text((nodes(edges_MMT(ii,1),1)+nodes(edges_MMT(ii,2),1))/2,(nodes(edges_MMT(ii,1),2)+nodes(edges_MMT(ii,2),2))/2,num2str(edges_MMT(ii,4)),'Fontsize',15);
            
        end
        for ii=1:nb.dirichlets
            line([nodes(dirichlets(ii,1),1) nodes(dirichlets(ii,2),1)],[nodes(dirichlets(ii,1),2) nodes(dirichlets(ii,2),2)]);
            text((nodes(dirichlets(ii,1),1)+nodes(dirichlets(ii,2),1))/2,(nodes(dirichlets(ii,1),2)+nodes(dirichlets(ii,2),2))/2,num2str(dirichlets(ii,4)),'Fontsize',15);
            
        end
        axis equal
        
    case 2
        
        
        
        
        figure
        hold on
        
        for ii=1:nb.loads
            line([nodes(loads(ii,1),1) nodes(loads(ii,2),1)],[nodes(loads(ii,1),2) nodes(loads(ii,2),2)],'Color','red');
        end
        for ii=1:nb.periodicity
            line([nodes(periodicity(ii,1),1) nodes(periodicity(ii,2),1)],[nodes(periodicity(ii,1),2) nodes(periodicity(ii,2),2)],'Color','k');
            
        end
        for ii=1:nb.MMT
            line([nodes(edges_MMT(ii,1),1) nodes(edges_MMT(ii,2),1)],[nodes(edges_MMT(ii,1),2) nodes(edges_MMT(ii,2),2)]);
            
        end
        for ii=1:nb.interfaces
            line([nodes(interfaces(ii,1),1) nodes(interfaces(ii,2),1)],[nodes(interfaces(ii,1),2) nodes(interfaces(ii,2),2)]);
            
        end
        
        for ii=1:nb.dirichlets
            line([nodes(dirichlets(ii,1),1) nodes(dirichlets(ii,2),1)],[nodes(dirichlets(ii,1),2) nodes(dirichlets(ii,2),2)]);
            
        end
        axis equal
        
        
end

end
