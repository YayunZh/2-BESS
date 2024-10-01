function y = func_plot_task_lsv(CL,ENV,ExpRes)
     %% Plot the Utilization ratio -- Initial design Converters VS Traditional
    figure(11);
    for i = 1:size(ExpRes.bat_uratio_energy,1)
        pic(i) = plot((ENV.Fir_Conv.Num*ENV.Sweep.Conv.e_lim)/(CL.Stat.Bat_num*CL.Bat{1}.curlim*CL.Bat{1}.volt),...
            mean(reshape(ExpRes.bat_uratio_energy(i,:,:),ENV.Sweep.Stat.Conv,ENV.Fir_Conv.MC_trial),2),'d-','linewidth',2);
        hold on;
    end
    % xlabel('Sum of Converter energy ratio in Traditional DPP');
    % ylabel('Bat Utiliation');
    xlabel('Converter Energy Processed / Battery Capacity');
    ylabel('Battery Capacity Utilization');
    legend(pic,'5% Capacity Variation','10% Capacity Variation','15% Capacity Variation',...
        '20% Capacity Variation','25% Capacity Variance','50% Capacity Variance',...
        '75% Capacity Variance','100% Capacity Variance','Fontsize',10);
    %legend(pic_trad,'Trad Bat energy Var = 0.01','Trad Bat energy Var = 0.02','Trad Bat energy Var = 0.03','Trad Bat energy Var = 0.04','Fontsize',10);
    grid on;
    grid minor;
    
        %% Plot
    figure(12);
    for i = 1:size(ExpRes.bat_uratio_energy,1)
        pic(i) = plot((ENV.Fir_Conv.Num*ENV.Sweep.Conv.e_lim)/(CL.Stat.Bat_num*CL.Bat{1}.curlim*CL.Bat{1}.volt),...
            100-100*mean(reshape(ExpRes.total_energy_process(i,:,:)*0.15./ExpRes.max_output_energy(i,:,:),ENV.Sweep.Stat.Conv,ENV.Fir_Conv.MC_trial),2),'d-','linewidth',2);
        hold on;
    end

    % xlabel('Sum of Converter energy Rate in Traditional DPP');
    % ylabel('Bat Utiliation');
    xlabel('Converter Energy Processed / Battery Capacity');
    ylabel('System Energy Efficiency (Assume 85% Converter Efficiency)');
    legend(pic,'5% Capacity Variance','10% Capacity Variance','15% Capacity Variance',...
        '20% Capacity Variance','25% Capacity Variance','50% Capacity Variance',...
        '75% Capacity Variance','100% Capacity Variance','Fontsize',10);
%    legend(pic_trad,'Trad Bat energy Var = 0.01','Trad Bat energy Var = 0.02','Trad Bat energy Var = 0.03','Trad Bat energy Var = 0.04','Fontsize',10);
%    xlim([0 1.5]);
    grid on;
    grid minor;
end