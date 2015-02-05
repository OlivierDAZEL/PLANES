%%%%% Coordinates of the edge's vertices

coord_edge(1:2,1)=nodes(internal(ie,1),1:2)';
coord_edge(1:2,2)=nodes(internal(ie,2),1:2)';

a=coord_edge(:,1);
b=coord_edge(:,2);

h=norm(b-a);
n_ab=(b-a)/h;

%%%%% Elements on both sides of the edge

e_1=internal(ie,3);
e_2=internal(ie,4);
c_1=mean(nodes(elements(e_1,:),1:2))';
c_2=mean(nodes(elements(e_2,:),1:2))';
e_edge=min([e_1,e_2]);

%%%%% vector normal pointing towards \Omega_e'

centre_temp=mean(nodes(elements(e_edge,:),1:2))'; 
centre_edge=(a+b)/2;
n_centre=centre_temp-centre_edge;
ne=normal_edge(coord_edge);
if (n_centre'*ne<0)
    ne=-ne;
end
nx=ne(1);
ny=ne(2);


valid_edge_internal=0;
if (element_label(e_1)==element_label(e_2))
    if ((element_label(e_1)==0)|(floor(element_label(e_edge)/1000)==2)|(floor(element_label(e_edge)/1000)==3))
        %disp('Lancement internal_fluid')
        internal_fluid
        valid_edge_internal=1;
    elseif (floor(element_label(e_edge)/1000)==4)
        %disp('Lancement internal_PEM')
        internal_PEM
        valid_edge_internal=1;
    end
else
    if (sum(floor(element_label(e_1)/1000)==[0 2 3]))*(sum(floor(element_label(e_2)/1000)==[0 2 3]))
        %disp('Lancement fluid_fluid')
        fluid_fluid
        valid_edge_internal=1;
    end
    if (sum(floor(element_label(e_1)/1000)==[0 2 3]))*(sum(floor(element_label(e_2)/1000)==[4]))
        %disp('Lancement fluid_PEM')
        fluid_PEM
        valid_edge_internal=1;
    end
    if (sum(floor(element_label(e_2)/1000)==[0 2 3]))*(sum(floor(element_label(e_1)/1000)==[4]))
        %disp('Lancement PEM_fluid')
        PEM_fluid
        valid_edge_internal=1;
    end

end

if valid_edge_internal==0
    disp('Stop in edge internal not a valid internal edge')
    jkljkkjklj
end


