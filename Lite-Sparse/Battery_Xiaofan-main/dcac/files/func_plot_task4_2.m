function y = func_plot_task4_2(CL,ENV,ExpRes, ENV_rand, ExpRes_rand)
%power_traj = 1.8*[0 0 0 0.2 0.5 1 1.1 1.2 1.33 1 0.8 0.67 0.8 1 1.1 1.2 1.1 0.95 0.8 0.6 0.4 0.33 0.33 0.33];
power_traj = 9*[0 0 0 0.2 0.5 1 1.1 1.2 1.33 1 0.8 0.67 0.8 1 1.1 1.2 1.1 0.95 0.8 0.6 0.4 0.33 0.33 0.33]/1.33;
scaling_factor = 80/9;
tl_ind = 1:24;
for k2 = 1:ENV.Var_Conv.MC_trial     % k_2 loop is for MC simu
    for i = 1:length(tl_ind)
        if (ExpRes.max_output_energy(3,2,k2) < power_traj(tl_ind(i)))
            ev_demand_gap_good(k2,tl_ind(i)) =  (power_traj(tl_ind(i)) - ExpRes.max_output_energy(3,2,k2))*scaling_factor;
        else
            ev_demand_gap_good(k2,tl_ind(i)) = 0;
        end
    end
end 

for k2 = 1:ENV_rand.Var_Conv.MC_trial     % k_2 loop is for MC simu
    for i = 1:length(tl_ind)
        if (ExpRes_rand.max_output_energy(3,2,k2) < power_traj(tl_ind(i)))
            ev_demand_gap_rand(k2,tl_ind(i)) =  (power_traj(tl_ind(i)) - ExpRes_rand.max_output_energy(3,2,k2))*scaling_factor;
        else
            ev_demand_gap_rand(k2,tl_ind(i)) = 0;
        end
    end
end

figure(41);
%for i = 1:size(ExpRes.bat_uratio_energy,1)
    pic(1) = plot(tl_ind,mean(ev_demand_gap_good,1),'d-','linewidth',2);
    hold on;
    pic(2) = plot(tl_ind,mean(ev_demand_gap_rand,1),'d-','linewidth',2);
    hold on;
    pic(3) = plot(tl_ind,power_traj*scaling_factor,'d-','linewidth',2);
    % xlabel('Sum of Converter Power Rate in Traditional DPP');
    % ylabel('Bat Utiliation');
    xlabel('Time (h)');
    ylabel('EV Charging Power (kW)');
    legend(pic,'Demand Gap in Optimized Power Flow','Demand Gap in Non-Optimized Power Flow','Demand','Fontsize',10);
%    legend(pic_trad,'Trad Bat Power Var = 0.01','Trad Bat Power Var = 0.02','Trad Bat Power Var = 0.03','Trad Bat Power Var = 0.04','Fontsize',10);
%    xlim([0 1.5]);
    grid on;
    grid minor;
end