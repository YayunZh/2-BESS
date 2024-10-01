function y = func_plot_task6(CL,ENV,TstBat)
%% Data Preparations
    TstBat_uratio_reorg = func_organize_4d(TstBat.bat_uratio_energy);
    TstBat_energy_process_reorg = func_organize_4d(TstBat.total_energy_process);
    TstBat_max_output_reorg = func_organize_4d(TstBat.max_output_energy);
    TstBat_energy_rating_reorg = func_organize_5d(TstBat.bat_energy_rating);
    
    for i = 1: ENV.Sweep.Stat.BatPopVar
        for j = 1:ENV.BatPopVar.MC_trial
            Bat_uratio_reorg(i,j) = TstBat.full_info{i,j}.bat_uratio;
        end
        Bat_uratio_reorg_mean(i) = mean(Bat_uratio_reorg(i,:));
    end    
   
    for i = 1: ENV.Sweep.Stat.BatPopVar
        for j = 1:ENV.BatPopVar.MC_trial
            Bat_processed_energy_reorg(i,j) = TstBat.full_info{i,j}.DesignRes.Total_energy_process;
        end
        Bat_processed_energy_reorg_mean = mean(Bat_processed_energy_reorg,2);
    end   

%     for i = 1:CL.Stat.Bat_num
%         InitialBat_energy(i) = CL.Bat{i}.volt*2;
%     end
%% Check the batery statistics 
    figure(60);
    Bat_statistics = func_bat_stat(TstBat.bat_energy_rating);
    subplot(3,3,1);
    histogram(Bat_statistics(1,1,:));
    title('Pop Var = 5%, Test Var = 5%');
    subplot(3,3,2);
    histogram(Bat_statistics(2,1,:));
    title('Pop Var = 10%, Test Var = 5%');
    subplot(3,3,3);
    histogram(Bat_statistics(3,1,:));
    title('Pop Var = 20%, Test Var = 5%');
    subplot(3,3,4);
    histogram(Bat_statistics(1,2,:));
    title('Pop Var = 5%, Test Var = 10%');
    subplot(3,3,5);
    histogram(Bat_statistics(2,2,:));
    title('Pop Var = 10%, Test Var = 10%');
    subplot(3,3,6);
    histogram(Bat_statistics(3,2,:));
    title('Pop Var = 20%, Test Var = 10%');
    subplot(3,3,7);
    histogram(Bat_statistics(1,3,:));
    title('Pop Var = 5%, Test Var = 20%');
    subplot(3,3,8);
    histogram(Bat_statistics(2,3,:));
    title('Pop Var = 10%, Test Var = 20%');
    subplot(3,3,9);
    histogram(Bat_statistics(3,3,:));
    title('Pop Var = 20%, Test Var = 20%');
    a = 1;
% Plot the Utilization ratio -- Initial design Converters VS Traditional
    figure(61);
    for i = 1:ENV.Sweep.Stat.BatPopVar 
        pic(i) = plot(100*ENV.Sweep.Bat{1}.qlim_test_var./ENV.Sweep.Bat{1}.qlim_mu,...
            100*Bat_uratio_reorg_mean(i) - 100*TstBat_uratio_reorg(i,:),'d-','linewidth',2);
        hold on;
    end
    xlabel('Sum of Converter Power ratio in Traditional DPP');
    ylabel('Bat Utilization');
    xlabel('Testing Uncertainty (%)');
    ylabel('Battery Energy Utilization Gap(%)');
    legend(pic,'5% Capacity Var','10% Capacity Var','20% Capacity Var','Fontsize',10);
    grid on;
    grid minor;
    
    figure(62);
   for i = 1:ENV.Sweep.Stat.BatPopVar 
        pic(i) = plot(100*ENV.Sweep.Bat{1}.qlim_test_var./ENV.Sweep.Bat{1}.qlim_mu,...
            100*TstBat_uratio_reorg(i,:),'d-','linewidth',2);
        hold on;
    end
    % xlabel('Sum of Converter Power ratio in Traditional DPP');
    % ylabel('Bat Utiliation');
    xlabel('Testing Uncertainty (%)');
    ylabel('Battery Energy Utilization(%)');
    legend(pic,'5% Capacity Var','10% Capacity Var','20% Capacity Var','Fontsize',10);
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
    function y = func_organize_4d(TstBat2)
        for ic1 = 1:size(TstBat2,1)
  %          for ic2 = 1:size(TstBat2,2)
                for ic3 = 1:size(TstBat2,3)
  %                  for ic4 = 1:size(TstBat2,4)
                         temp = reshape(TstBat2(ic1,:,ic3,:),1,size(TstBat2,2)*size(TstBat2,4));
                         TstBat_reorg2(ic1,ic3) = mean(temp);
  %                  end
                end
  %          end
        end
        y = TstBat_reorg2;
    end 

function y = func_bat_stat(TstBat2)
    for ic1 = 1:size(TstBat2,1)
%          for ic2 = 1:size(TstBat2,2)
            for ic3 = 1:size(TstBat2,3)
%                  for ic4 = 1:size(TstBat2,4)
                 temp(ic1,ic3,:) = reshape(TstBat2(ic1,:,ic3,:,:),1,size(TstBat2,2)*size(TstBat2,4)*size(TstBat2,5));
%                  end
            end
%          end
    end
    y = temp;
end 

    function y = func_organize_5d(TstBat2)
        for ic1 = 1:size(TstBat2,1)
  %          for ic2 = 1:size(TstBat2,2)
                for ic3 = 1:size(TstBat2,3)
  %                  for ic4 = 1:size(TstBat2,4)
                        for ic5 = 1:size(TstBat2,5)
                             temp(ic5,:) = reshape(TstBat2(ic1,:,ic3,:,ic5),1,size(TstBat2,2)*size(TstBat2,4));
                        end
                        TstBat_reorg2(ic1,ic3) = mean(sum(temp,1));
  %                  end
                end
  %          end
        end
        y = TstBat_reorg2;
    end 
end