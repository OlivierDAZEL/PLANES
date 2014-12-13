


fid=fopen(nom_fichier_msh,'r');
nb_noeuds=fscanf(fid,'%i',1);
nb_elements=fscanf(fid,'%i',1);
nb_edges=fscanf(fid,'%i',1);

for ii=1:nb_noeuds
    vcor(ii,1)=fscanf(fid,'%f',1);
    vcor(ii,2)=fscanf(fid,'%f',1);
    vcor(ii,3)=0;
    node_label(ii)=fscanf(fid,'%f',1);
end
for ii=1:nb_elements
    kconec(ii,1)=fscanf(fid,'%i',1);
    kconec(ii,2)=fscanf(fid,'%i',1);
    kconec(ii,3)=fscanf(fid,'%i',1);
    element_label(ii)=fscanf(fid,'%f',1);
end
for ii=1:nb_edges
    edge(ii,1)=fscanf(fid,'%i',1);
    edge(ii,2)=fscanf(fid,'%i',1);
    edge(ii,3)=fscanf(fid,'%i',1);
end
fclose(fid);




% figure
% 
% 
%  for ii=1:nb_edges
%      line([vcor(edge(ii,1),1) vcor(edge(ii,2),1)],[vcor(edge(ii,1),2) vcor(edge(ii,2),2)],'Color','g');
%      text(mean(vcor(edge(ii,1:2),1)),mean(vcor(edge(ii,1:2),2)),num2str(edge(ii,3)),'Fontsize',15);
%  end
% 
% axis equal









% vcor=load('swap/vcor.txt');
% kconec=load('swap/kconec.txt');
% edge=load('swap/edge.txt');
% element_label=load('swap/element_label.txt');

nb_elements=size(kconec,1);
nb_noeuds=size(vcor,1);
nb_edges=size(edge,1);
figure
hold on

for ii=1:nb_elements
    line([vcor(kconec(ii,1),1) vcor(kconec(ii,2),1)],[vcor(kconec(ii,1),2) vcor(kconec(ii,2),2)],'Color','r');
    line([vcor(kconec(ii,2),1) vcor(kconec(ii,3),1)],[vcor(kconec(ii,2),2) vcor(kconec(ii,3),2)],'Color','r');
    line([vcor(kconec(ii,3),1) vcor(kconec(ii,4),1)],[vcor(kconec(ii,3),2) vcor(kconec(ii,4),2)],'Color','r');
    line([vcor(kconec(ii,4),1) vcor(kconec(ii,5),1)],[vcor(kconec(ii,4),2) vcor(kconec(ii,5),2)],'Color','r');
    line([vcor(kconec(ii,5),1) vcor(kconec(ii,6),1)],[vcor(kconec(ii,5),2) vcor(kconec(ii,6),2)],'Color','r');
    line([vcor(kconec(ii,6),1) vcor(kconec(ii,1),1)],[vcor(kconec(ii,6),2) vcor(kconec(ii,1),2)],'Color','r');    
    %text(mean(vcor(kconec(ii,:),1)),mean(vcor(kconec(ii,:),2)),num2str(ii),'Fontsize',15);
end

plot(vcor(:,1),vcor(:,2),'.','Markersize',20);
for ii=1:nb_noeuds
    text(vcor(ii,1),vcor(ii,2),num2str(ii),'Fontsize',15);
end

figure
hold on

for ii=1:nb_elements
    line([vcor(kconec(ii,1),1) vcor(kconec(ii,2),1)],[vcor(kconec(ii,1),2) vcor(kconec(ii,2),2)],'Color','r');
    line([vcor(kconec(ii,2),1) vcor(kconec(ii,3),1)],[vcor(kconec(ii,2),2) vcor(kconec(ii,3),2)],'Color','r');
    line([vcor(kconec(ii,3),1) vcor(kconec(ii,4),1)],[vcor(kconec(ii,3),2) vcor(kconec(ii,4),2)],'Color','r');
    line([vcor(kconec(ii,4),1) vcor(kconec(ii,5),1)],[vcor(kconec(ii,4),2) vcor(kconec(ii,5),2)],'Color','r');
    line([vcor(kconec(ii,5),1) vcor(kconec(ii,6),1)],[vcor(kconec(ii,5),2) vcor(kconec(ii,6),2)],'Color','r');
    line([vcor(kconec(ii,6),1) vcor(kconec(ii,1),1)],[vcor(kconec(ii,6),2) vcor(kconec(ii,1),2)],'Color','r');   
    text(mean(vcor(kconec(ii,:),1)),mean(vcor(kconec(ii,:),2)),num2str(element_label(ii)),'Fontsize',15);
end
axis equal
figure

% plot(vcor(:,1),vcor(:,2),'.','Markersize',20);
 for ii=1:nb_edges
%     if edge(ii,3)==3
     line([vcor(edge(ii,1),1) vcor(edge(ii,6),1)],[vcor(edge(ii,1),2) vcor(edge(ii,6),2)],'Color','g');
     line([vcor(edge(ii,6),1) vcor(edge(ii,2),1)],[vcor(edge(ii,6),2) vcor(edge(ii,2),2)],'Color','m');     
     text(mean(vcor(edge(ii,1:2),1)),mean(vcor(edge(ii,1:2),2)),num2str(edge(ii,3)),'Fontsize',15);
 %    end
 end

axis equal