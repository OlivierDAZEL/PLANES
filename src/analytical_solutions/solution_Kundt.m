



X_an=1/sin(k_air*ly);


x_air=linspace(-ly,0,1000);


sigma_a=air.K*k_air*X_an*cos(k_air*x_air);


if plot_profiles==1
    
    
    % figure (1001)
    % hold on
    
    
    
    figure(10002)
    hold on
    plot(x_air+ly,abs(sigma_a))
    
    % figure(1000)
    % subplot(3,3,1)
    % plot(x_EF,abs(v_s),'r')
    % subplot(3,3,3)
    % hold on
    % plot(x_EF,abs(v_t),'r')
    % plot(x_air,abs(v_x_air))
    % subplot(3,3,5)
    % plot(x_EF,abs(sigma_xx+sigma_yy)/2,'r')
    % subplot(3,3,8)
    % hold on
    % plot(x_EF,abs(sigma_p),'r')
    % plot(x_air,abs(sigma_a))
    % subplot(3,3,2)
    % plot(x_EF,abs(0),'r')
    % subplot(3,3,4)
    % hold on
    % plot(x_air,abs(0))
    % plot(x_EF,abs(0),'r')
    % subplot(3,3,6)
    % plot(x_EF,abs(0),'r')
    % subplot(3,3,7)
    % plot(x_EF,abs(sigma_xx-sigma_yy)/2,'r')
    %
    % %
    % %
    % %
    % % % figure(1001)
    % % % subplot(3,3,1)
    % % % plot(x_EF+d_poreux+d_air,angle(v_s),'r')
    % % % subplot(3,3,3)
    % % % hold on
    % % % plot(x_EF+d_poreux+d_air,angle(v_t),'r')
    % % % plot(x_air+d_poreux+d_air,angle(v_x_air))
    % % % subplot(3,3,5)
    % % % plot(x_EF+d_poreux+d_air,angle(sigma_xx+sigma_yy),'r')
    % % % subplot(3,3,8)
    % % % hold on
    % % % plot(x_EF+d_poreux+d_air,angle(sigma_p),'r')
    % % % plot(x_air+d_poreux+d_air,angle(sigma_a))
    % % % subplot(3,3,2)
    % % % plot(x_EF+d_poreux+d_air,angle(0),'r')
    % % % subplot(3,3,4)
    % % % hold on
    % % % plot(x_air+d_poreux+d_air,angle(0))
    % % % plot(x_EF+d_poreux+d_air,angle(0),'r')
    % % % subplot(3,3,6)
    % % % plot(x_EF+d_poreux+d_air,angle(0),'r')
    % % % subplot(3,3,7)
    % % % plot(x_EF+d_poreux+d_air,angle(sigma_xx-sigma_yy),'r')
    % %
end

% 
% 
% for e_edge=1:nb_elements
%     centre_elem=mean(nodes(elements(e_edge,:),1:2))';
%     parameter_element
%     e_plus=exp(+1j*k_e*(centre_elem(2)-(d_air+d_EF)));
%     abs(e_plus)
%     e_moins=1/e_plus;
%     
%     if floor(element_label(e_edge)/1000)==2
%         
%         
%         
%         q_theorique(dof_start_element(e_edge)+1)=e_plus*omega*A/2;
%         q_theorique(dof_start_element(e_edge)  )=e_moins*omega*A/2;
%         
%         
%     end
%     if  floor(element_label(e_edge)/1000)==0
%         
% 
%         q_theorique(dof_start_element(e_edge)+1)=e_plus*omega*((C/2-D/(2*i)));
%         q_theorique(dof_start_element(e_edge)  )=e_moins*omega*((C/2+D/(2*i)));
%         
%     end
% end
% 
% 
% figure (21)
% hold on
% 
% 
% plot(abs(X),'r.')
% plot(abs(q_theorique),'k+')
% %plot(abs(q_theorique*X(27)/q_theorique(27)),'k+')
% % 
% % figure (22)
% % hold on
% % plot(angle(X),'r.')
% % plot(angle(-q_theorique),'k+')
% 
