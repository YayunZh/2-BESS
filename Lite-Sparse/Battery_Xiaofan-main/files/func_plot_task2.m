function y = func_plot_task2(CL,ENV,TstBat)
%% Data Preparations
    TstBat_uratio_reorg = func_organize_4d(TstBat.bat_uratio_power);
    TstBat_power_process_reorg = func_organize_4d(TstBat.total_power_process);
    TstBat_max_output_reorg = func_organize_4d(TstBat.max_output_power);
    TstBat_power_rating_reorg = func_organize_5d(TstBat.bat_power_rating);
    
    for i = 1: ENV.Sweep.Stat.BatPopVar
        for j = 1:ENV.BatPopVar.MC_trial
            Bat_uratio_reorg(i,j) = TstBat.full_info{i,j}.bat_uratio;
        end
        Bat_uratio_reorg_mean(i) = mean(Bat_uratio_reorg(i,:));
    end    
   
    for i = 1: ENV.Sweep.Stat.BatPopVar
        for j = 1:ENV.BatPopVar.MC_trial
            Bat_processed_power_reorg(i,j) = TstBat.full_info{i,j}.DesignRes.Total_power_process;
        end
        Bat_processed_power_reorg_mean = mean(Bat_processed_power_reorg,2);
    end   
    % Plot the Utilization Rate -- Initial design Converters
%    figure(21)
%     pic1 = plot(ENV.Sweep.Bat{1}.var(end),...
%          mean(reshape(ExpRes.urate(end,ceil(ENV.Sweep.Stat.Conv*0.1),:),1,ENV.Var_Conv.MC_trial),2),'d-','linewidth',2);
%     hold on;
    % Plot the Utilization Rate -- Redesign Converters
%     pic2 = plot(ENV2.Sweep.Bat{1}.var,...
%             mean(reshape(ExpRes_redesign.urate(:,1,:),ENV2.Sweep.Stat.Bat,ENV2.Var_Conv.MC_trial),2),'--','linewidth',2);
%    for i = 1:size(ExpRes_redesign.urate,1)-1
%         pic2(i) = plot(1-ENV2.Sweep.TstBatT,...
%             mean(reshape(ExpRes_redesign.max_power(i,:,:)./ExpRes_redesign.urate(i,:,:),ENV2.Sweep.Stat.TstBatT, ENV2.Tst_Bat.MC_trial),2),'d-','linewidth',2);
%            pic2(i) = plot(1-ENV2.Sweep.TstBatT,...
%               mean(reshape(ExpRes_redesign.urate(i,:,:),ENV2.Sweep.Stat.TstBatT, ENV2.Tst_Bat.MC_trial),2),'d-','linewidth',2);   
%     hold on;
%     end
%    xlabel('Testing Time');
%     ylabel('Max Output Power');
%    ylabel('Battery Utilization');
%    legend([pic1, pic2],'Original Design','Redesign','Fontsize',10);
    %legend( pic2,'Var Pop = 0.05','Var Pop = 0.15','Var Pop = 0.25','Var Pop = 0.35','Fontsize',10);
%    legend( pic2,'Var = 0.39','Var = 0.37','Var = 0.3','Fontsize',10);
%     legend( pic2,'Var = 0.39','Var = 0.37','Var = 0.3','Fontsize',10);
%     grid on;
%     grid minor;

  %% Plot the Utilization ratio -- Initial design Converters VS Traditional
%     figure(21);
%     for i = 1:ENV.Sweep.Stat.BatPopVar 
%         pic(i) = plot(100*ENV.Sweep.Bat{1}.curlim_test_var./ENV.Sweep.Bat{1}.curlim_mu,...
%             100*TstBat.full_info{i}.bat_uratio - 100*mean(reshape(TstBat.bat_uratio_power(i,:,:),...
%             size(TstBat.bat_uratio_power(i,:,:),2), size(TstBat.bat_uratio_power(i,:,:),3)),2),'d-','linewidth',2);
%         hold on;
%     end
%     % xlabel('Sum of Converter Power ratio in Traditional DPP');
%     % ylabel('Bat Utiliation');
%     xlabel('Testing Uncertainty (%)');
%     ylabel('Battery Utilization Gap(%)');
%     % legend(pic,'Bat power Var = 5% Averaged power','Bat power Var = 10% Averaged power','Bat power Var = 20% Averaged power','Fontsize',10);
%     grid on;
%     grid minor;
% 
%     figure(22);
%     bat_uratio_mat = 100*reshape(TstBat.bat_uratio_power(2,:,:), size(TstBat.bat_uratio_power(i,:,:),2),  size(TstBat.bat_uratio_power(i,:,:),3));
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
%% Data Preparations
    TstBat_uratio_reorg = func_organize_4d(TstBat.bat_uratio_power);
    TstBat_power_process_reorg = func_organize_4d(TstBat.total_power_process);
    TstBat_max_output_reorg = func_organize_4d(TstBat.max_output_power);
    TstBat_power_rating_reorg = func_organize_5d(TstBat.bat_power_rating);
    
    for i = 1: ENV.Sweep.Stat.BatPopVar
        for j = 1:ENV.BatPopVar.MC_trial
            Bat_uratio_reorg(i,j) = TstBat.full_info{i,j}.bat_uratio;
        end
        Bat_uratio_reorg_mean = mean(Bat_uratio_reorg,2);
    end    
   
    for i = 1: ENV.Sweep.Stat.BatPopVar
        for j = 1:ENV.BatPopVar.MC_trial
            Bat_processed_power_reorg(i,j) = TstBat.full_info{i,j}.DesignRes.Total_power_process;
        end
        Bat_processed_power_reorg_mean = mean(Bat_processed_power_reorg,2);
    end   

%     for i = 1:CL.Stat.Bat_num
%         InitialBat_power(i) = CL.Bat{i}.volt*2;
%     end
    %% Check the batery statistics 
    figure(60);
    Bat_statistics = func_bat_stat(TstBat.bat_power_rating);
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
    
    figure(61);
    for i = 1:ENV.Sweep.Stat.BatPopVar 
        pic(i) = plot(100*ENV.Sweep.Bat{1}.curlim_test_var./ENV.Sweep.Bat{1}.curlim_mu,...
            100*Bat_uratio_reorg_mean(i) - 100*TstBat_uratio_reorg(i,:),'d-','linewidth',2);
        hold on;
    end
    % xlabel('Sum of Converter Power ratio in Traditional DPP');
    % ylabel('Bat Utiliation');
    xlabel('Testing Uncertainty (%)');
    ylabel('Battery Power Utilization Gap(%)');
    legend(pic,'5% Power Var','10% Power Var','20% Power Var','Fontsize',10);
    grid on;
    grid minor;
    
    figure(62);
    for i = 1:ENV.Sweep.Stat.BatPopVar 
        pic(i) = plot(100*ENV.Sweep.Bat{1}.curlim_test_var./ENV.Sweep.Bat{1}.curlim_mu,...
            100*TstBat_uratio_reorg(i,:),'d-','linewidth',2);
        hold on;
    end
    % xlabel('Sum of Converter Power ratio in Traditional DPP');
    % ylabel('Bat Utiliation');
    xlabel('Testing Uncertainty (%)');
    ylabel('Battery Energy Utilization(%)');
    legend(pic,'5% Power Var','10% Power Var','20% Power Var','Fontsize',10);
    grid on;
    grid minor;

%     figure(62);
%     bat_uratio_mat = 100*reshape(TstBat.bat_uratio_power(2,:,:),  size(TstBat.bat_uratio_power(i,:,:),2),  size(TstBat.bat_uratio_power(i,:,:),3));
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
