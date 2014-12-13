isvalidof=zeros(3*nb_nodes,1);

for ie=1:nb_elements

	typ=floor(element_label(ie)/1000);
	switch typ
		case {0,2,3,8} %! Acoustic/EF/limep
					isvalidof(3*(elements(ie,1:6)-1)+3)=1;
        case {1} %! Elastic solid
					isvalidof(3*(elements(ie,1:6)-1)+1)=1;
					isvalidof(3*(elements(ie,1:6)-1)+2)=1;

         case {4,5}	%! PEM
					isvalidof(3*(elements(ie,1:6)-1)+1)=1;
					isvalidof(3*(elements(ie,1:6)-1)+2)=1;
					isvalidof(3*(elements(ie,1:6)-1)+3)=1;
        otherwise
				disp('Attention element sans nature connue')
                stop
    end
end   



for ie=1:nb_dirichlets
	if (dirichlets(ie,4)==5) % Sliding
		xx=abs(nodes(dirichlets(ie,1),1)-nodes(dirichlets(ie,2),1));
		yy=abs(nodes(dirichlets(ie,1),2)-nodes(dirichlets(ie,2),2));
		if (xx>yy) 
			isvalidof(3*(dirichlets(ie,1)-1)+2)=0;
			isvalidof(3*(dirichlets(ie,2)-1)+2)=0;
			isvalidof(3*(dirichlets(ie,6)-1)+2)=0;
        else
			isvalidof(3*(dirichlets(ie,1)-1)+1)=0;
			isvalidof(3*(dirichlets(ie,2)-1)+1)=0;
			isvalidof(3*(dirichlets(ie,6)-1)+1)=0;
        end
    end
	if (dirichlets(ie,4)==6) % Bonded
        isvalidof(3*(dirichlets(ie,1)-1)+1)=0;
		isvalidof(3*(dirichlets(ie,1)-1)+2)=0;
		isvalidof(3*(dirichlets(ie,2)-1)+1)=0;
		isvalidof(3*(dirichlets(ie,2)-1)+2)=0;
		isvalidof(3*(dirichlets(ie,6)-1)+1)=0;
		isvalidof(3*(dirichlets(ie,6)-1)+2)=0;
    end
	if (dirichlets(ie,4)==1) 
		isvalidof(3*(dirichlets(ie,1)-1)+1)=0;
		isvalidof(3*(dirichlets(ie,1)-1)+2)=0;
		isvalidof(3*(dirichlets(ie,2)-1)+1)=0;
		isvalidof(3*(dirichlets(ie,2)-1)+2)=0;
		isvalidof(3*(dirichlets(ie,6)-1)+1)=0;
		isvalidof(3*(dirichlets(ie,6)-1)+2)=0;
    end	
end





dof_A=zeros(3*nb_nodes,1);


dof_back=[];



itemp=0;
for ii=1:3*nb_nodes
	if (isvalidof(ii)) 
		itemp=itemp+1;
		dof_A(ii)=itemp;
		dof_back(itemp)=ii;
    end
end

nb_dof_FEM=itemp;
list_dof_valid=find(isvalidof);

