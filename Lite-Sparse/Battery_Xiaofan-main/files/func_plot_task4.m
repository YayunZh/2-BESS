function y = func_plot_task4(CL,ENV,ExpRes)
    %% Plot the Utilization Rate -- Initial design Converters VS Traditional
    figure(41);
    for i = 1:size(ExpRes.bat_uratio_energy,1)
        pic(i) = plot(ENV.Var_Conv.Num*ENV.Sweep.Conv.e_lim./(sum(ENV.Avg_Conv.e_lim_vec)),...
            mean(reshape(ExpRes.bat_uratio_energy(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2),'d-','linewidth',2);
        hold on;
    end
    % xlabel('Sum of Converter Power Rate in Traditional DPP');
    % ylabel('Bat Utiliation');
    xlabel('\lambda (Converter Scaling)');
    ylabel('Bat Energy Capacity Utilization');
    legend(pic,'5% Capacity Variance','10% Capacity Variance','15% Capacity Variance',...
        '20% Capacity Variance','Fontsize',10);
%    legend(pic_trad,'Trad Bat Power Var = 0.01','Trad Bat Power Var = 0.02','Trad Bat Power Var = 0.03','Trad Bat Power Var = 0.04','Fontsize',10);
%    xlim([0 1.5]);
    grid on;
    grid minor;
    
        %% Plot
    figure(42);
    for i = 1:size(ExpRes.bat_uratio_energy,1)
        pic(i) = plot(ENV.Var_Conv.Num*ENV.Sweep.Conv.e_lim./(sum(ENV.Avg_Conv.e_lim_vec)),...
            mean(reshape(ExpRes.total_energy_process(i,:,:)./ExpRes.max_output_energy(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2),'d-','linewidth',2);
        hold on;
    end


    % xlabel('Sum of Converter Power Rate in Traditional DPP');
    % ylabel('Bat Utiliation');
    xlabel('\lambda (Converter Scaling)');
    ylabel('Processed Energy / Output Energy');
    legend(pic,'5% Capacity Variance','10% Capacity Variance','15% Capacity Variance',...
        '20% Capacity Variance','Fontsize',10);
%    legend(pic_trad,'Trad Bat Power Var = 0.01','Trad Bat Power Var = 0.02','Trad Bat Power Var = 0.03','Trad Bat Power Var = 0.04','Fontsize',10);
%    xlim([0 1.5]);
    grid on;
    grid minor;
    
    %% Plot
    figure(43);
    for i = 1:size(ExpRes.bat_uratio_energy,1)
        pic(i) = plot(ENV.Var_Conv.Num*ENV.Sweep.Conv.e_lim./(sum(ENV.Avg_Conv.e_lim_vec)),...
            1-mean(reshape(ExpRes.total_energy_process(i,:,:)*0.85./ExpRes.max_output_energy(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2),'d-','linewidth',2);
        hold on;
    end
    % xlabel('Sum of Converter Power Rate in Traditional DPP');
    % ylabel('Bat Utiliation');
    xlabel('\lambda (Converter Scaling)');
    ylabel('System Energy Efficiency ((Assume 85% Converter Efficiency)');
    legend(pic,'5% Capacity Variance','10% Capacity Variance','15% Capacity Variance',...
        '20% Capacity Variance','Fontsize',10);
%    legend(pic_trad,'Trad Bat Power Var = 0.01','Trad Bat Power Var = 0.02','Trad Bat Power Var = 0.03','Trad Bat Power Var = 0.04','Fontsize',10);
%    xlim([0 1.5]);
    grid on;
    grid minor;
end