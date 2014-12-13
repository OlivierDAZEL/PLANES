nb_noeuds_lin=size(vcor,1);
nb_segments=3*size(kconec,1);
nb_noeuds=size(vcor,1);
nb_elements=size(kconec,1);
nb_edge=size(edge,1);

segments=zeros(2,nb_segments);
segment_ok=ones(nb_segments,1);

for ii=1:nb_elements
    segments(1,1+3*(ii-1))=kconec(ii,1);
    segments(2,1+3*(ii-1))=kconec(ii,2);
    segments(1,2+3*(ii-1))=kconec(ii,2);
    segments(2,2+3*(ii-1))=kconec(ii,3);
    segments(1,3+3*(ii-1))=kconec(ii,3);
    segments(2,3+3*(ii-1))=kconec(ii,1);
end



for ii=1:nb_segments
    if (segments(1,ii)>segments(2,ii))
        jj=segments(2,ii);
        segments(2,ii)=segments(1,ii);
        segments(1,ii)=jj;
    end
end


for ii=1:(nb_segments-1)
    if (segment_ok(ii)==1)
        for jj=ii+1:nb_segments
            if (segments(1,ii)==segments(1,jj))&(segments(2,ii)==segments(2,jj))
                segment_ok(jj)=0;
            end
        end
    end
end



nvcbis=0;
for ii=1:nb_segments
    if (segment_ok(ii)==1)
        nvcbis=nvcbis+1;
    end
end


vcl=zeros(nb_noeuds+nvcbis,3);
til=zeros(nb_elements,6);


for ii=1:nb_elements
    for jj=1:3
        til(ii,jj)=kconec(ii,jj);
    end
end
for ii=1:nb_noeuds
    vcl(ii,1)=vcor(ii,1);
    vcl(ii,2)=vcor(ii,2);
end

for ii=1:nb_elements
    til(ii,5)=til(ii,3);
    til(ii,3)=til(ii,2);
    til(ii,2)=0;
end


nvc_add=0;
for ii=1:nb_elements
    jj=1;
    noeud_1=til(ii,1);
    noeud_2=til(ii,3);

    if (segment_ok(jj+3*(ii-1))==1)
        nvc_add=nvc_add+1;
        vcl(nb_noeuds+nvc_add,1)=(vcl(noeud_1,1)+vcl(noeud_2,1))/2;
        vcl(nb_noeuds+nvc_add,2)=(vcl(noeud_1,2)+vcl(noeud_2,2))/2;
        til(ii,2)=nb_noeuds+nvc_add;
    else
        x_temp=(vcl(noeud_1,1)+vcl(noeud_2,1))/2;
        y_temp=(vcl(noeud_1,2)+vcl(noeud_2,2))/2;
        [min_noeud i_noeud]=min(abs((vcl(:,1)-x_temp).^2+(vcl(:,2)-y_temp).^2));
        til(ii,2)=i_noeud;
    end
    jj=2;
    noeud_1=til(ii,3);
    noeud_2=til(ii,5);

    if (segment_ok(jj+3*(ii-1))==1)
        nvc_add=nvc_add+1;
        vcl(nb_noeuds+nvc_add,1)=(vcl(noeud_1,1)+vcl(noeud_2,1))/2;
        vcl(nb_noeuds+nvc_add,2)=(vcl(noeud_1,2)+vcl(noeud_2,2))/2;
        til(ii,4)=nb_noeuds+nvc_add;
    else
        x_temp=(vcl(noeud_1,1)+vcl(noeud_2,1))/2;
        y_temp=(vcl(noeud_1,2)+vcl(noeud_2,2))/2;
        [min_noeud i_noeud]=min(abs((vcl(:,1)-x_temp).^2+(vcl(:,2)-y_temp).^2));
        til(ii,4)=i_noeud;
    end
    jj=3;
    noeud_1=til(ii,5);
    noeud_2=til(ii,1);

    if (segment_ok(jj+3*(ii-1))==1)
        nvc_add=nvc_add+1;
        vcl(nb_noeuds+nvc_add,1)=(vcl(noeud_1,1)+vcl(noeud_2,1))/2;
        vcl(nb_noeuds+nvc_add,2)=(vcl(noeud_1,2)+vcl(noeud_2,2))/2;
        til(ii,6)=nb_noeuds+nvc_add;
    else
        x_temp=(vcl(noeud_1,1)+vcl(noeud_2,1))/2;
        y_temp=(vcl(noeud_1,2)+vcl(noeud_2,2))/2;
        [min_noeud i_noeud]=min(abs((vcl(:,1)-x_temp).^2+(vcl(:,2)-y_temp).^2));
        til(ii,6)=i_noeud;
    end


end


nb_nodes=nb_noeuds+nvcbis;
kconec=til;
vcor=vcl;
edge2=edge;

for ie=1:nb_edge
    noeud_1=edge(ie,1);
    noeud_2=edge(ie,2);
    %On cherche le noeud milieu
    x_temp=(vcor(noeud_1,1)+vcor(noeud_2,1))/2;
    y_temp=(vcor(noeud_1,2)+vcor(noeud_2,2))/2;
    [min_noeud i_noeud]=min(abs((vcor(:,1)-x_temp).^2+(vcor(:,2)-y_temp).^2));
    edge2(ie,6)=i_noeud;
end


nb_edges=size(edge2,1);
edge=edge2;


clear x_temp
clear y_temp
clear i_noeud
clear min_noeud
clear nvc_add
clear nvcbis
clear nb_segments
clear segments_ok
clear segments
