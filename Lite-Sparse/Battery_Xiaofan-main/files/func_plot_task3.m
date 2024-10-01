function y = func_plot_task3(ExpRes,CL,ENV)
%     figure(31);
%     for i = 1:size(ExpRes.urate,1)
%         pic(i) = plot(ENV.Var_Conv.Num*ENV.Sweep.Conv.p_lim + sum(ENV.Avg_Conv.p_lim),...
%             mean(reshape(ExpRes.bat_eff(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2),'d-','linewidth',2);
%         hold on;
%     end
%     % xlabel('Sum of Converter Power Rate in Traditional DPP');
%     % ylabel('Bat Utiliation');
%     xlabel('Sum of Avg and Var Converter Power Rate');
%     ylabel('Battery Effciency');
%     legend(pic,'Bat Res Var = 50% Mean','Bat Res Var = 34% Mean','Bat Res Var = 16% Mean','Bat Res Var = 0% Mean','Fontsize',10);
%     grid on;
%     grid minor;
% 
%     figure(32);
%     for i = 1:size(ExpRes.urate,1)
%         pic(i) = plot(ENV.Var_Conv.Num*ENV.Sweep.Conv.p_lim + sum(ENV.Avg_Conv.p_lim),...
%             mean(reshape(ExpRes.bat_eff(i,:,:).*ExpRes.bat_eff_error(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2),'d-','linewidth',2);
%         hold on;
%     end
%     % xlabel('Sum of Converter Power Rate in Traditional DPP');
%     % ylabel('Bat Utiliation');
%     xlabel('Sum of Avg and Var Converter Power Rate');
%     ylabel('Battery Effciency Error');
%     legend(pic,'Bat Res Var = 50% Mean','Bat Res Var = 34% Mean','Bat Res Var = 16% Mean','Bat Res Var = 0% Mean','Fontsize',10);
%     grid on;
%     grid minor;

 figure(51);
    for i = 1:ENV.Sweep.Stat.Loss
        pic(i) = plot((ENV.Var_Conv.Num*ENV.Sweep.Conv.p_lim+sum(ENV.Avg_Conv.p_lim_vec))/9,...
            100*mean(reshape(ExpRes.bat_power_eff_error(i,:,:),ENV.Sweep.Stat.Conv,ENV.PostEff.MC_trial),2),'d-','linewidth',2);
        hold on;
    end
    % xlabel('Sum of Converter Power Rate in Traditional DPP');
    % ylabel('Bat Utiliation');
    % ylim([0 2]);
    title('BESS Efficiency with 95% Efficiency Batteries')
    xlabel('Normalized Aggregate Converter Rating');
    ylabel('BESS Power Efficiency Error (%)');
    legend(pic,'Bat Eff Std = 5% Averaged Bat Eff','Bat Eff Std = 10% Averaged Bat Eff','Bat Eff Std = 15% Averaged Bat Eff','Bat Eff Std = 20% Averaged Bat Eff','Fontsize',10);
    grid on;
    grid minor;
    
  figure(52);
    for i = 1:ENV.Sweep.Stat.Loss
        pic(i) = plot((ENV.Var_Conv.Num*ENV.Sweep.Conv.p_lim+sum(ENV.Avg_Conv.p_lim_vec))/9,...
            100*mean(reshape(ExpRes.bat_power_eff(i,:,:),ENV.Sweep.Stat.Conv,ENV.PostEff.MC_trial),2),'d-','linewidth',2);
        hold on;
    end
    % xlabel('Sum of Converter Power Rate in Traditional DPP');
    % ylabel('Bat Utiliation');
    ylim([90 100]);
    title('BESS Efficiency with 95% Efficiency Batteries');
    xlabel('Normalized Aggregate Converter Rating');
    ylabel('BESS Power Efficiency (%)');
    legend(pic,'Bat Eff Std = 5% Averaged Bat Eff','Bat Eff Std = 10% Averaged Bat Eff','Bat Eff Std = 15% Averaged Bat Eff','Bat Eff Std = 20% Averaged Bat Eff','Fontsize',10);
    grid on;
    grid minor;
    
    %%
%     figure(52);
%     for i = 1:size(ExpRes.bat_uratio_power,1)
%         pic(i) = plot(ENV.Var_Conv.Num*ENV.Sweep.Conv.p_lim./(sum(ENV.Avg_Conv.p_lim)),...
%             mean(reshape(ExpRes.bat_power_eff(i,:,:),ENV.Sweep.Stat.Conv,ENV.PostEff.MC_trial),2),'d-','linewidth',2);
%         hold on;
%     end
%     % xlabel('Sum of Converter Power Rate in Traditional DPP');
%     % ylabel('Bat Utiliation');
%  %   ylim([0 2]);
%     xlabel('\lambda (Converter scaling)');
%     ylabel('Battery Power Efficiency (%)');
%     legend(pic,'Bat Eff Std = 2.5% Averaged Bat Eff','Bat Eff Std = 5% Averaged Bat Eff','Bat Eff Std = 7.5% Averaged Bat Eff','Bat Eff Std = 10% Averaged Bat Eff','Fontsize',10);
%     grid on;
%     grid minor;

