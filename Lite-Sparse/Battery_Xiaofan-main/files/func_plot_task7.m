function y = func_plot_task7(CL,ENV,DiagnosisBat)
    %% Plot the Utilization ratio -- Initial design Converters VS Traditional
%% Data Preparations
    DiagnosisBat_uratio_reorg = func_organize_4d(DiagnosisBat.bat_uratio_energy);
    DiagnosisBat_energy_process_reorg = func_organize_4d(DiagnosisBat.total_energy_process);
    DiagnosisBat_max_output_reorg = func_organize_4d(DiagnosisBat.max_output_energy);
    DiagnosisBat_energy_rating_reorg = func_organize_5d(DiagnosisBat.bat_energy_rating);
    
    for i = 1: ENV.Sweep.Stat.BatPopVar
        for j = 1:ENV.BatPopVar.MC_trial
            InitialBat_uratio_reorg(i,j) = DiagnosisBat.full_info{i,j}.bat_uratio;
        end
    end    
   
    for i = 1: ENV.Sweep.Stat.BatPopVar
        for j = 1:ENV.BatPopVar.MC_trial
            InitialBat_processed_energy_reorg(i,j) = DiagnosisBat.full_info{i,j}.DesignRes.Total_energy_process;
        end
    end   

    for i = 1:CL.Stat.Bat_num
        InitialBat_energy(i) = CL.Bat{i}.volt*2;
    end
    
%% Plotting, show to Anna
% for multiple x and y axis, please refer to 
% https://www.mathworks.com/help/matlab/creating_plots/graph-with-multiple-x-axes-and-y-axes.html
%     figure(71);
% %     for i = 2:ENV.Sweep.Stat.BatPopVar 
% %    for i = 1:1
% %         pic(i) = plot([0 100 300 500],...
% %             [0,100*mean(InitialBat_uratio_reorg(i,:)) - 100*DiagnosisBat_uratio_reorg(i,:)],... 
% %                 'd-','linewidth',2);
% %         hold on;
% %     end
%     pic(1) = plot([0 100 300 500],...
%             [0,100*mean(InitialBat_uratio_reorg(2,:)) - 100*DiagnosisBat_uratio_reorg(2,:)],... 
%                 'd-','linewidth',2,'color','r');
%             hold on;
%     pic(2) = plot([0 100 300 500],...
%             [0,100*mean(InitialBat_uratio_reorg(2,:)) - 100*DiagnosisBat_uratio_reorg(2,:)],... 
%                 'd--','linewidth',2,'color','r');
%             hold on;
%     pic(3) = plot([0 100 300 500],...
%             [0,100*mean(InitialBat_uratio_reorg(3,:)) - 100*DiagnosisBat_uratio_reorg(3,:)],... 
%                 'd-.','linewidth',2,'color','r');
%     xlabel('Cycle Numbers');
%     ylabel('Battery Utilization Gap (%)');
%     %legend(pic,'5% Energy Variation','10% Energy Variation','20% Energy Variation','Fontsize',10);
%     legend(pic,'5% Energy Variation','10% Energy Variation','20% Energy Variation','Fontsize',10);
%      ylim([-2,12]);
% %     grid on;
% %     grid minor;
%     ax1 = gca; % current axes
%     ax1.XColor = 'r';
%     ax1.YColor = 'r';
%     ax1_pos = ax1.Position; % position of first axes
%     ax2 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','None');
%      
%     for i = 1:1 
%         pic2(i) = line(100*[0,ENV.Sweep.Bat{1}.qlim_diagnosis_var./(ENV.Sweep.Bat{1}.qlim_mu)'],...
%            100*[1,DiagnosisBat_energy_rating_reorg(i,:)/sum(InitialBat_energy)],... 
%                 'Parent',ax2,'Color','k','LineStyle','--','LineWidth',3);
%     end
%     xlabel('Degradation Variation (%)');
%     ylabel('Average Battery Capacity (%)');
%  %   legend(pic,'Bat Energy Var = 5% Averaged Energy','Bat Energy Var = 10% Averaged Energy','Bat Energy Var = 20% Averaged Energy','Fontsize',10);
%     grid on;
%     grid minor;
    
    
 %% Plotting, show to myself
    figure(71);   
%     for i = 1:ENV.Sweep.Stat.BatPopVar 
%         pic(i) = plot([0,100*ENV.Sweep.Bat{1}.qlim_diagnosis_var./(ENV.Sweep.Bat{1}.qlim_mu)'],...
%            [mean(InitialBat_processed_energy_reorg(i,:)),DiagnosisBat_energy_process_reorg(i,:)]/sum(InitialBat_energy),... 
%                 '--','linewidth',2);
%         hold on;
%     end  
     for i = 1:ENV.Sweep.Stat.BatPopVar 
        pic(i) = plot([0,100*ENV.Sweep.Bat{1}.qlim_diagnosis_var./(ENV.Sweep.Bat{1}.qlim_mu)'],...
             [0,100*mean(InitialBat_uratio_reorg(i,:)) - 100*DiagnosisBat_uratio_reorg(i,:)],... 
                'd-','linewidth',2);
        hold on;
    end
    xlabel('Cycle Number');
    ylabel('Processed Energy/Initial Battery Energy');
    legend(pic,'5% Energy Variation','10% Energy Variation','20% Energy Variation','Fontsize',10);
    grid on;
    grid minor;
    
%     figure(62);
%     bat_uratio_mat = 100*reshape(TstBat.bat_uratio_energy(2,:,:),  size(TstBat.bat_uratio_energy(i,:,:),2),  size(TstBat.bat_uratio_energy(i,:,:),3));
%     subplot(1,3,1);
%     histogram(bat_uratio_mat(1,:),30,'Normalization','probability');
%     title('5% Testing Uncertainty');
%     xlabel('Battery Utilization (%)');
%     ylabel('Probability');
%     grid on;
%     grid minor;
%     subplot(1,3,2);
%     histogram(bat_uratio_mat(2,:),30,'Normalization','probability');
%     title('10% Testing Uncertainty');
%     xlabel('Battery Utilization (%)');
%     ylabel('Probability');
%     grid on;
%     grid minor;
%     subplot(1,3,3);
%     histogram(bat_uratio_mat(3,:),30,'Normalization','probability');
%     title('20% Testing Uncertainty');
%     xlabel('Battery Utilization (%)');
%     ylabel('Probability');
%     grid on;
%     grid minor;
    

%% User Functions
    function y = func_organize_4d(DiagnosisBat2)
        for ic1 = 1:size(DiagnosisBat2,1)
  %          for ic2 = 1:size(DiagnosisBat2,2)
                for ic3 = 1:size(DiagnosisBat2,3)
  %                  for ic4 = 1:size(DiagnosisBat2,4)
                         temp = reshape(DiagnosisBat2(ic1,:,ic3,:),1,size(DiagnosisBat2,2)*size(DiagnosisBat2,4));
                         DiagnosisBat_reorg2(ic1,ic3) = mean(temp);
  %                  end
                end
  %          end
        end
        y = DiagnosisBat_reorg2;
    end 

    function y = func_organize_5d(DiagnosisBat2)
        for ic1 = 1:size(DiagnosisBat2,1)
  %          for ic2 = 1:size(DiagnosisBat2,2)
                for ic3 = 1:size(DiagnosisBat2,3)
  %                  for ic4 = 1:size(DiagnosisBat2,4)
                        for ic5 = 1:size(DiagnosisBat2,5)
                             temp(ic5,:) = reshape(DiagnosisBat2(ic1,:,ic3,:,ic5),1,size(DiagnosisBat2,2)*size(DiagnosisBat2,4));
                        end
                        DiagnosisBat_reorg2(ic1,ic3) = mean(sum(temp,1));
  %                  end
                end
  %          end
        end
        y = DiagnosisBat_reorg2;
    end 
    
end