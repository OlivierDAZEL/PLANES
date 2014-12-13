





for ie=1:nb_elements
    if floor(element_label(ie)/1000)==0
        ondes_element(ie)=1;   % Air
    elseif floor(element_label(ie)/1000)==1
        ondes_element(ie)=2;   % Milieu solide elastique
    elseif floor(element_label(ie)/1000)==2
        ondes_element(ie)=1;   % Fluide equivalent
     elseif floor(element_label(ie)/1000)==3
        ondes_element(ie)=1;   % Limp
      elseif floor(element_label(ie)/1000)==4
        ondes_element(ie)=3;   % Fluide equivalent
      elseif floor(element_label(ie)/1000)==8
        ondes_element(ie)=1;   % Fluide equivalent  
    else
        disp('Subroutine import_mesh')
        disp('Unknwon type of element')
        break
    end
end

angle_element=nb_theta*ones(nb_elements,1);

dof_start_element(1)=1; %
for ie=2:nb_elements
   dof_start_element(ie)=dof_start_element(ie-1)+ondes_element(ie-1)*angle_element(ie-1);
end

