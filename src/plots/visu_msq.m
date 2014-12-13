% visu_msq.m
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