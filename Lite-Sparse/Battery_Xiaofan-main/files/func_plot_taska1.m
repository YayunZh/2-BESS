function y = func_plot_taska1(CL,ENV, ExpRes)
    %% Plot the Utilization ratio -- Initial design Converters VS Traditional
    figure(11);
    pic = plot((ENV.Var_Conv.Num*ENV.Sweep.Conv.p_lim + sum(ENV.Avg_Conv.p_lim))/(CL.Stat.Bat_num*CL.Bat{1}.curlim*CL.Bat{1}.volt),...
        mean(ExpRes.bat_uratio_power(:,:),2),'d-','linewidth',2);
    % xlabel('Sum of Converter Power ratio in Traditional DPP');
    % ylabel('Bat Utiliation');
    xlabel('Converter Power Rating / Battery Power Rating');
    ylabel('Bat Power Capability Utilization');
  %  legend(pic,'5% Power Capability Variation','10% Power Capability Variation','15% Power Capability Variation','20% Power Capability Variation','Fontsize',10);
    %legend(pic_trad,'Trad Bat Power Var = 0.01','Trad Bat Power Var = 0.02','Trad Bat Power Var = 0.03','Trad Bat Power Var = 0.04','Fontsize',10);
    grid on;
    grid minor;

    %% Plot the Battery Utilization ratio and Converter Utilization ratio Tradeoff
%     figure(12);
%     for i = 1:size(ExpRes.bat_uratio_power,1)
%         pic(i) = plot(mean(reshape(ExpRes.uratio_conv(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2),...
%             mean(reshape(ExpRes.bat_uratio_power(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2),'d-','linewidth',2);
%         hold on;
%     end
%     xlabel('Converter Utiliation');
%     ylabel('Bat Utiliation');
%     legend(pic,'Bat Power Var = 0.01','Bat Power Var = 0.02','Bat Power Var = 0.03','Bat Power Var = 0.04','Fontsize',10);
%     grid on;
%     grid minor;
 %% Plot the Processed Power
%     figure(13);
%     for i = 1:size(ExpRes.process_power,1)
%         pic(i) = plot(ENV.Sweep.Conv.p_lim, mean(reshape(ExpRes.process_power(i,:,:), ENV.Sweep.Stat.Conv, ENV.Var_Conv.MC_trial),2),'s-','linewidth',2);
%         hold on;
%     end
%     xlabel('Var Converter Power');
%     ylabel('Processed Power');
%     legend(pic,'Bat Power Var = 0.01','Bat Power Var = 0.02','Bat Power Var = 0.03','Bat Power Var = 0.04','Fontsize',10);
%     grid on;
%     grid minor;
% 
  %% Plot the Maximum Output Power
%     figure(14);
%     for i = 1:size(ExpRes.max_power,1)
%         pic(i) = plot(ENV.Sweep.Conv.p_lim, mean(reshape(ExpRes.max_power(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2),'*-','linewidth',2);
%         hold on;
%     end
%     xlabel('Var Converter Power');
%     ylabel('Max Output Power');
%     legend(pic,'Bat Power Var = 0.01','Bat Power Var = 0.02','Bat Power Var = 0.03','Bat Power Var = 0.04','Fontsize',10);
%     grid on;
%     grid minor;
% 
  %% Plot the partial power ratio
%     figure(15);
%     ppp_ratio = cell(size(ExpRes.max_power,1),size(ExpRes.max_power,2));
%     ppp_ratio_mean = [];
%     for k0 = 1:size(ExpRes.max_power,1)
%         for k1 = 1:size(ExpRes.max_power,2)
%             for k2 = 1:size(ExpRes.max_power,3)
%                 if(ExpRes.max_power(k0,k1,k2) ~= 0)
%                     ppp_ratio{k0,k1}(end+1) = ExpRes.process_power(k0,k1,k2)./ExpRes.max_power(k0,k1,k2);
%                 else
%                     ppp_ratio{k0,k1}(end+1) = 1;   % indicating we have to use full power processing
%                 end
%                 ppp_ratio_mean(k0,k1) = mean(ppp_ratio{k0,k1});
%             end
%         end
%     end
%     for i = 1:size(ExpRes.max_power,1)
%         pic(i) = plot(ENV.Sweep.Conv.p_lim, ppp_ratio_mean(i,:),'d--','linewidth',2);
%         hold on;
%     end
%     xlabel('Var Converter Power');
%     ylabel('Partial Power Ratio');
%     legend(pic,'Bat Power Var = 0.01','Bat Power Var = 0.02','Bat Power Var = 0.03','Bat Power Var = 0.04','Fontsize',10);
%     grid on;
%     grid minor;
end