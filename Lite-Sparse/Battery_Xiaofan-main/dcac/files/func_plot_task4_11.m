function y = func_plot_task4_11(ENV,ExpTrajRes,ExpTrajRes_trad)

scalefac1 = 150/36;
% ExpTrajRes_typ = ExpTrajRes{1,1,1}{1,1};
% ExpTrajRes_typ_trad = ExpTrajRes_trad{1,1,1}{1,1};

% ExpTrajRes_typ.Pbessout = ExpectedTrajPbessout (ExpTrajRes, ENV);
% ExpTrajRes_typ_trad.Pbessout = ExpectedTrajPbessout (ExpTrajRes_trad, ENV);
for kn1 = 1:ENV.Sweep.Stat.EV_Mu                       % k_n1 loop for ev demand mean      
    for k0 = 1:ENV.Sweep.Stat.EV_Var                   % k_0 loop for ev demand variance      
        for k1 = 1:ENV.Sweep.Stat.Td_Mu                % k_1 loop for Tdelay mu 
            LP_Charging_time_matrix(kn1,k0,k1) = ExpectedLPChargingTime_method2(ExpTrajRes, ENV, kn1, k0, k1);
            LP_Charging_time_trad_matrix(kn1,k0,k1) = ExpectedLPChargingTime_method2(ExpTrajRes_trad, ENV, kn1, k0, k1);
        end
    end
end

figure();
for k1 = 1:ENV.Sweep.Stat.Td_Mu
    pic(k1) = plot(ENV.Sweep.EV_Var.e_var_normalized, 60*LP_Charging_time_matrix(2,:,k1),'-s','linewidth',2);
    hold on;
end
for k1 = 1:ENV.Sweep.Stat.Td_Mu
    pic2(k1) = plot(ENV.Sweep.EV_Var.e_var_normalized, 60*LP_Charging_time_matrix(1,:,k1),'--s','linewidth',2);
    hold on;
end
xlabel('EV Demand Deviation');
ylabel('Curtailed Charging Time (min)');
legend([pic,pic2],'$\bar{T}_d$ = 0.5 h $\bar{E}_{ev}$ = 50 kWh',...
'$\bar{T}_d$ = 1.5 h $\bar{E}_{ev}$ = 50 kWh',... 
'$\bar{T}_d$ = 2.5 h $\bar{E}_{ev}$ = 50 kWh',...
'$\bar{T}_d$ = 0.5 h $\bar{E}_{ev}$ = 33 kWh',...
'$\bar{T}_d$ = 1.5 h $\bar{E}_{ev}$ = 33 kWh',...
'$\bar{T}_d$ = 2.5 h $\bar{E}_{ev}$ = 33 kWh','Interpreter','latex');
title('LS-HiPPP');
ylim([0,0.7*60]);
grid on;
grid minor;

figure();
for k1 = 1:ENV.Sweep.Stat.Td_Mu
    pic(k1) = plot(ENV.Sweep.EV_Var.e_var_normalized, 60*LP_Charging_time_trad_matrix(2,:,k1),'-d','linewidth',2);
    hold on;
end
for k1 = 1:ENV.Sweep.Stat.Td_Mu
    pic2(k1) = plot(ENV.Sweep.EV_Var.e_var_normalized, 60*LP_Charging_time_trad_matrix(1,:,k1),'--d','linewidth',2);
    hold on;
end
xlabel('EV Demand Deviation');
ylabel('Curtailed Charging Time (min)');
legend([pic,pic2],'$\bar{T}_d$ = 0.5 h $\bar{E}_{ev}$ = 50 kWh',...
'$\bar{T}_d$ = 1.5 h $\bar{E}_{ev}$ = 50 kWh',... 
'$\bar{T}_d$ = 2.5 h $\bar{E}_{ev}$ = 50 kWh',...
'$\bar{T}_d$ = 0.5 h $\bar{E}_{ev}$ = 33 kWh',...
'$\bar{T}_d$ = 1.5 h $\bar{E}_{ev}$ = 33 kWh',...
'$\bar{T}_d$ = 2.5 h $\bar{E}_{ev}$ = 33 kWh','Interpreter','latex');
% legend([pic,pic2],'AvgTd = 0.5h AvgEV = 50 kWh','AvgTd = 1.5h AvgEV = 50 kWh','AvgTd = 2.5h AvgEV = 50 kWh',...
% 'AvgTd = 0.5h AvgEV = 33 kWh','AvgTd = 1.5h AvgEV = 33 kWh','AvgTd = 2.5h AvgEV = 33 kWh'...
% );
title('C-PPP');
ylim([0,0.7*60]);
grid on;
grid minor;

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

function lp_charging_time_mean = ExpectedLPChargingTime(ExpTrajRes, ENV, kn1, k0, k1)
    lp_charging_time_vec = [];
    for k2 = 1:ENV.Var_Conv.MC_trial
        for k3 = 1:ENV.Sweep.Stat.EV_MC_trial
            for k4 = 1:ENV.Sweep.Stat.Td_MC_trial
                LPcharging_time_cycle = [];
                for i_sweep = 1:length(ExpTrajRes{kn1,k0,k1,k2}{k3,k4}.Pevin.value)
                    if (ExpTrajRes{kn1,k0,k1,k2}{k3,k4}.Pevin.value(i_sweep) > 0) && (ExpTrajRes{kn1,k0,k1,k2}{k3,k4}.Pevin.value(i_sweep) < 36)
                        LPcharging_time_cycle = [LPcharging_time_cycle, ExpTrajRes{kn1,k0,k1,k2}{k3,k4}.Pevin.delta_time(i_sweep)];
                    end
                end
                if(isempty(LPcharging_time_cycle))
                    lp_charging_time_cycle_mean = 0;
                else
                    lp_charging_time_cycle_mean = mean(LPcharging_time_cycle);
                end
                lp_charging_time_vec = [lp_charging_time_vec,lp_charging_time_cycle_mean];
            end
        end
    end
    lp_charging_time_mean = mean(lp_charging_time_vec);
end

function lp_charging_time_mean = ExpectedLPChargingTime_method2(ExpTrajRes, ENV, kn1, k0, k1)
    lp_charging_time_vec = [];
    for k2 = 1:ENV.Var_Conv.MC_trial
        for k3 = 1:ENV.Sweep.Stat.EV_MC_trial
            for k4 = 1:ENV.Sweep.Stat.Td_MC_trial
                lp_charging_time_vec = [lp_charging_time_vec,ExpTrajRes{kn1,k0,k1,k2}{k3,k4}.CurtChargeTime.value];
            end
        end
    end
    lp_charging_time_mean = mean(lp_charging_time_vec);
end


end

