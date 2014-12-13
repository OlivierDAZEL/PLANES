
x=[];
y=[];
figure
hold on
for ie=1:nb_elements
    
    
    c=X_EF_fortran(kconec(ie,:));
    
    vertices=[vcor(kconec(ie,:),:)'];
    
    %changement de base (cart to triangle)
    coord=[];
    for i=1:length(vertices(1,:))
        coord=[coord vertices(1:2,i)-vertices(1:2,1)];
    end
    coeffksi=zeros(2,1);
    coeffeta=zeros(2,1);
    M=[coord(:,3)' 0 0; 0 0 coord(:,5)'; coord(:,5)' 0 0; 0 0 coord(:,3)'];
    coeffbase=inv(M)*[1;1;0;0];
    coeffksi=coeffbase(1:2);
    coeffeta=coeffbase(3:4);
    
    
    faces = [1 2 3]';
    vert=coord(:,[1 3 5]);
    for i=1:1
        [vert, faces]=linearSubdivision(vert, faces);
    end
    
    ksi=vert'*coeffksi;
    eta=vert'*coeffeta;
    lambda=1-eta-ksi;
    mat_ksi_eta=[-lambda.*(1-2*lambda) 4*ksi.*lambda -ksi.*(1-2*ksi) 4*ksi.*eta -eta.*(1-2*eta) 4*eta.*lambda];
    
    val_int=mat_ksi_eta*c;
    
    for i_faces=1:size(faces,2)
        figure(0006)
        hold on
        plot(vert(1,faces(:,i_faces))+vertices(1,1)',abs(val_int(faces(:,i_faces))),'.');
    end
    
    
    
    
    
end
title('Pression EF (module)')




