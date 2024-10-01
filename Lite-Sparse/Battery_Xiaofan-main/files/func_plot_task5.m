function y = func_plot_task5(ExpRes,CL,ENV)
    figure(51);
    for i = 1:size(ExpRes.bat_uratio_energy,1)
        pic(i) = plot(ENV.Var_Conv.Num*ENV.Sweep.Conv.e_lim./(sum(ENV.Avg_Conv.e_lim)),...
            100*mean(reshape(ExpRes.bat_energy_eff_error(i,:,:),ENV.Sweep.Stat.Conv,ENV.PostEff.MC_trial),2),'d-','linewidth',2);
        hold on;
    end
    % xlabel('Sum of Converter Power Rate in Traditional DPP');
    % ylabel('Bat Utiliation');
 %   ylim([0 2]);
    xlabel('\lambda (Converter scaling)');
    ylabel('System Energy Effciency Error (%)');
    legend(pic,'Bat Eff Std = 2.5% Averaged Bat Eff','Bat Eff Std = 5% Averaged Bat Eff','Bat Eff Std = 7.5% Averaged Bat Eff','Bat Eff Std = 10% Averaged Bat Eff','Fontsize',10);
    grid on;
    grid minor;
    
    %%
    figure(52);
    for i = 1:size(ExpRes.bat_uratio_energy,1)
        pic(i) = plot(ENV.Var_Conv.Num*ENV.Sweep.Conv.e_lim./(sum(ENV.Avg_Conv.e_lim)),...
            mean(reshape(ExpRes.bat_energy_eff(i,:,:),ENV.Sweep.Stat.Conv,ENV.PostEff.MC_trial),2),'d-','linewidth',2);
        hold on;
    end
    % xlabel('Sum of Converter Power Rate in Traditional DPP');
    % ylabel('Bat Utiliation');
 %   ylim([0 2]);
    xlabel('\lambda (Converter scaling)');
    ylabel('Error on efficiency (%)');
    legend(pic,'Bat Eff Std = 2.5% Averaged Bat Eff','Bat Eff Std = 5% Averaged Bat Eff','Bat Eff Std = 7.5% Averaged Bat Eff','Bat Eff Std = 10% Averaged Bat Eff','Fontsize',10);
    grid on;
    grid minor;