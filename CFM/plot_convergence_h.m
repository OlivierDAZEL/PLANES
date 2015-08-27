clear all
close all


load('../Projects/CFM_square/CFM_square.error')



     figure;
     loglog(CFM_square(:,1),CFM_square(:,5),'-')
    
   
    
    
%          figure;
%     semilogx(nb_dof_list,NNZ_list(:,:,nwave)'./nb_dof_list,'-')
% 
%     hold on
%     semilogx(FEM(:,3),FEM(:,5)./FEM(:,3),'.','Markersize',30)
%     set(gca,'LineWidth',1,'FontSize',13);
%     set(get(gca,'Children'),'LineWidth',1.5);
%     xlabel('number of dof');
%     legend('N_w=4','N_w=6','N_w=8','N_w=10','N_w=12','N_w=14','N_w=16','FEM','Location','NorthEast');
%     ylabel('NNz./nb_dof');
%     print('-depsc2', sprintf('figures/error_p_h_%s_rapports.eps', wave_name{nwave}));
    
    
    
%     figure;
%     loglog(1./H_list,Er_L2_v_list(:,:,nwave),'-')
%     axis([1 30 1.e-8 30])
%     set(gca,'LineWidth',1,'FontSize',13);
%     set(get(gca,'Children'),'LineWidth',1.5);
%     xlabel('1/h');
%     ylabel('L2 error on p');
%     print('-depsc2', sprintf('figures/error_v_h_%s.eps', wave_name{nwave}));
    
%     figure;
%     loglog(1./H_list,Er_L2_vs_list(:,:,nwave),'-')
%     axis([1 30 1.e-8 30])
%     set(gca,'LineWidth',1,'FontSize',13);
%     set(get(gca,'Children'),'LineWidth',1.5);
%     xlabel('1/h');
%     ylabel('L2 error on p');
%     print('-depsc2', sprintf('figures/error_vs_h_%s.eps', wave_name{nwave}));    




