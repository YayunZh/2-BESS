function y = func_plot_task4_9(ENV,ExpTrajRes,ExpTrajRes_trad)

scalefac1 = 150/36;
% ExpTrajRes_typ = ExpTrajRes{1,1,1}{1,1};
% ExpTrajRes_typ_trad = ExpTrajRes_trad{1,1,1}{1,1};

% ExpTrajRes_typ.Pbessout = ExpectedTrajPbessout (ExpTrajRes, ENV);
% ExpTrajRes_typ_trad.Pbessout = ExpectedTrajPbessout (ExpTrajRes_trad, ENV);
for kn1 = 1:ENV.Sweep.Stat.EV_Mu                       % k_n1 loop for ev demand mean      
    for k0 = 1:ENV.Sweep.Stat.EV_Var                   % k_0 loop for ev demand variance      
        for k1 = 1:ENV.Sweep.Stat.Td_Mu                % k_1 loop for Tdelay mu 
            Charging_time_matrix(kn1,k0,k1) = ExpectedChargingTime(ExpTrajRes, ENV, kn1, k0, k1);
            Charging_time_trad_matrix(kn1,k0,k1) = ExpectedChargingTime(ExpTrajRes_trad, ENV, kn1, k0, k1);
            Individual_Charging_time_matrix(kn1,k0,k1) = ExpectedTrajIndividualChargingTime(ExpTrajRes, ENV, kn1, k0, k1);
            Individual_Charging_time_trad_matrix(kn1,k0,k1) = ExpectedTrajIndividualChargingTime(ExpTrajRes_trad, ENV, kn1, k0, k1);
        end
    end
end

%% Total Charging Time
% figure();
% for k0 = 1:ENV.Sweep.Stat.EV_Var
%     pic(k0) = plot(ENV.Sweep.Td_Mu.mu, reshape(Charging_time_matrix(2,k0,:),1,ENV.Sweep.Stat.Td_Mu),'-s','linewidth',2);
%     hold on;
% end
% for k0 = 1:ENV.Sweep.Stat.EV_Var
%     pic2(k0) = plot(ENV.Sweep.Td_Mu.mu, reshape(Charging_time_matrix(1,k0,:),1,ENV.Sweep.Stat.Td_Mu),'--s','linewidth',2);
%     hold on;
% end
% xlabel('Td Mean (h)');
% ylabel('Charging Time (h)');
% legend([pic,pic2],'StdEV = 10% AvgEV = 50 kWh','StdEV = 30% AvgEV = 50 kWh','StdEV = 50% AvgEV = 50 kWh',...
% 'StdEV = 10% AvgEV = 33 kWh','StdEV = 30% AvgEV = 33 kWh','StdEV = 50% AvgEV = 33 kWh'...
% );
% title('LS-HiPPP');
% %ylim([0.6,1]);
% grid on;
% grid minor;
% 
% figure();
% for k0 = 1:ENV.Sweep.Stat.EV_Var
%     pic(k0) = plot(ENV.Sweep.Td_Mu.mu, reshape(Charging_time_trad_matrix(2,k0,:),1,ENV.Sweep.Stat.Td_Mu),'-s','linewidth',2);
%     hold on;
% end
% for k0 = 1:ENV.Sweep.Stat.EV_Var
%     pic2(k0) = plot(ENV.Sweep.Td_Mu.mu, reshape(Charging_time_trad_matrix(1,k0,:),1,ENV.Sweep.Stat.Td_Mu),'--s','linewidth',2);
%     hold on;
% end
% xlabel('Td Mean (h)');
% ylabel('Charging Time (h)');
% legend([pic,pic2],'StdEV = 10% AvgEV = 50 kWh','StdEV = 30% AvgEV = 50 kWh','StdEV = 50% AvgEV = 50 kWh',...
% 'StdEV = 10% AvgEV = 33 kWh','StdEV = 30% AvgEV = 33 kWh','StdEV = 50% AvgEV = 33 kWh'...
% );
% title('PPP');
% %ylim([0.6,1]);
% grid on;
% grid minor;

%% Individual Charging Time

figure();
for k1 = 1:ENV.Sweep.Stat.Td_Mu
    pic(k1) = plot(ENV.Sweep.EV_Var.e_var_normalized, 60*Individual_Charging_time_matrix(2,:,k1),'-s','linewidth',2);
    hold on;
end
for k1 = 1:ENV.Sweep.Stat.Td_Mu
    pic2(k1) = plot(ENV.Sweep.EV_Var.e_var_normalized, 60*Individual_Charging_time_matrix(1,:,k1),'--s','linewidth',2);
    hold on;
end
xlabel('E_{CG} Normalized Energy Gap');
ylabel('Individual Charging Time (min)');
legend([pic,pic2],'AvgTd = 0.5h AvgEV = 50 kWh','AvgTd = 1.5h AvgEV = 50 kWh','AvgTd = 2.5h AvgEV = 50 kWh',...
'AvgTd = 0.5h AvgEV = 33 kWh','AvgTd = 1.5h AvgEV = 33 kWh','AvgTd = 2.5h AvgEV = 33 kWh'...
);
title('LS-HiPPP');
ylim([0,0.7*60]);
grid on;
grid minor;

figure();
for k1 = 1:ENV.Sweep.Stat.Td_Mu
    pic(k1) = plot(ENV.Sweep.EV_Var.e_var_normalized, 60*Individual_Charging_time_trad_matrix(2,:,k1),'-d','linewidth',2);
    hold on;
end
for k1 = 1:ENV.Sweep.Stat.Td_Mu
    pic2(k1) = plot(ENV.Sweep.EV_Var.e_var_normalized, 60*Individual_Charging_time_trad_matrix(1,:,k1),'--d','linewidth',2);
    hold on;
end
xlabel('E_{CG} Normalized Energy Gap');
ylabel('Individual Charging Time (min)');
legend([pic,pic2],'AvgTd = 0.5h AvgEV = 50 kWh','AvgTd = 1.5h AvgEV = 50 kWh','AvgTd = 2.5h AvgEV = 50 kWh',...
'AvgTd = 0.5h AvgEV = 33 kWh','AvgTd = 1.5h AvgEV = 33 kWh','AvgTd = 2.5h AvgEV = 33 kWh'...
);
title('C-PPP');
ylim([0,0.7*60]);
grid on;
grid minor;


%% User Defined Functions
function uratio_energy_mean = ExpectedTrajuratio_energy(ExpTrajRes, ENV, kn1, k0, k1)
    uratio_energy_vec = [];
    for k2 = 1:ENV.Var_Conv.MC_trial
        for k3 = 1:ENV.Sweep.Stat.EV_MC_trial
            for k4 = 1:ENV.Sweep.Stat.Td_MC_trial
                uratio_energy_vec = [uratio_energy_vec, ExpTrajRes{kn1,k0,k1,k2}{k3,k4}.uratio_energy.value(2:end)];
            end
        end
    end
    uratio_energy_mean = mean(uratio_energy_vec);
end

function avg_chargingtime_mean = ExpectedTrajIndividualChargingTime(ExpTrajRes, ENV, kn1, k0, k1)
    avg_chargingtime_vec = [];
    for k2 = 1:ENV.Var_Conv.MC_trial
        for k3 = 1:ENV.Sweep.Stat.EV_MC_trial
            for k4 = 1:ENV.Sweep.Stat.Td_MC_trial
                avg_chargingtime_vec = [avg_chargingtime_vec, ExpTrajRes{kn1,k0,k1,k2}{k3,k4}.TotalChargeTime.value(2:end)];
            end
        end
    end
    avg_chargingtime_mean = mean(avg_chargingtime_vec);
end

function charging_time_mean = ExpectedChargingTime(ExpTrajRes, ENV, kn1, k0, k1)
    charging_time_vec = [];
    for k2 = 1:ENV.Var_Conv.MC_trial
        for k3 = 1:ENV.Sweep.Stat.EV_MC_trial
            for k4 = 1:ENV.Sweep.Stat.Td_MC_trial
                charging_time = 0;
                for i_sweep = 1:length(ExpTrajRes{kn1,k0,k1,k2}{k3,k4}.Pevin.value)
                    if (ExpTrajRes{kn1,k0,k1,k2}{k3,k4}.Pevin.value(i_sweep) > 0)
                        charging_time = charging_time + ExpTrajRes{kn1,k0,k1,k2}{k3,k4}.Pevin.delta_time(i_sweep);
                    end
                end
                charging_time_vec = [charging_time_vec, charging_time];
            end
        end
    end
    charging_time_mean = mean(charging_time_vec);
end
% function ExpTrajRes_Pbessout = ExpectedTrajPbessout(ExpTrajRes, ENV)
%     
%     Pbessout_matrix = [];
%     for i = 1:ENV.Var_Conv.MC_trial
%         Pbessout_matrix = [Pbessout_matrix; ExpTrajRes{1,1,i}{1,1}.Pbessout.value];
%     end
% 
%     ExpTrajRes_Pbessout.delta_time = ExpTrajRes{1,1,1}{1,1}.Pbessout.delta_time;
%     ExpTrajRes_Pbessout.abs_time = ExpTrajRes{1,1,1}{1,1}.Pbessout.abs_time;
%     ExpTrajRes_Pbessout.upper = max(Pbessout_matrix,[],1);
%     ExpTrajRes_Pbessout.lower = min(Pbessout_matrix,[],1);
%     ExpTrajRes_Pbessout.mean = mean(Pbessout_matrix,1);
%     
% 
% end


end

