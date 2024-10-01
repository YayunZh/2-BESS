function y = func_plot_task4_8(ENV,ExpTrajRes,ExpTrajRes_trad)

scalefac1 = 150/36;
% ExpTrajRes_typ = ExpTrajRes{1,1,1}{1,1};
% ExpTrajRes_typ_trad = ExpTrajRes_trad{1,1,1}{1,1};

% ExpTrajRes_typ.Pbessout = ExpectedTrajPbessout (ExpTrajRes, ENV);
% ExpTrajRes_typ_trad.Pbessout = ExpectedTrajPbessout (ExpTrajRes_trad, ENV);
for kn1 = 1:ENV.Sweep.Stat.EV_Mu                       % k_n1 loop for ev demand mean      
    for k0 = 1:ENV.Sweep.Stat.EV_Var                   % k_0 loop for ev demand variance      
        for k1 = 1:ENV.Sweep.Stat.Td_Mu                % k_1 loop for Tdelay mu 
            Charge_energy_matrix(kn1,k0,k1) = ExpectedChargerEnergy(ExpTrajRes, ENV, kn1, k0, k1);
            Charge_energy_trad_matrix(kn1,k0,k1) = ExpectedChargerEnergy(ExpTrajRes_trad, ENV, kn1, k0, k1);
        end
    end
end

figure();
for k0 = 1:ENV.Sweep.Stat.EV_Var
    pic(k0) = plot(ENV.Sweep.Td_Mu.mu, scalefac1*reshape(Charge_energy_matrix(2,k0,:),1,ENV.Sweep.Stat.Td_Mu),'-s','linewidth',2);
    hold on;
end
for k0 = 1:ENV.Sweep.Stat.EV_Var
    pic2(k0) = plot(ENV.Sweep.Td_Mu.mu, scalefac1*reshape(Charge_energy_matrix(1,k0,:),1,ENV.Sweep.Stat.Td_Mu),'--s','linewidth',2);
    hold on;
end
xlabel('Td Mean (h)');
ylabel('Individual Charging Time  (kWh)');
legend([pic,pic2],'StdEV = 10% AvgEV = 50 kWh','StdEV = 30% AvgEV = 50 kWh','StdEV = 50% AvgEV = 50 kWh',...
'StdEV = 10% AvgEV = 33 kWh','StdEV = 30% AvgEV = 33 kWh','StdEV = 50% AvgEV = 33 kWh'...
);
title('LS-HiPPP');
%ylim([0.6,1]);
grid on;
grid minor;

figure();
for k0 = 1:ENV.Sweep.Stat.EV_Var
    pic(k0) = plot(ENV.Sweep.Td_Mu.mu, scalefac1*reshape(Charge_energy_trad_matrix(2,k0,:),1,ENV.Sweep.Stat.Td_Mu),'-s','linewidth',2);
    hold on;
end
for k0 = 1:ENV.Sweep.Stat.EV_Var
    pic2(k0) = plot(ENV.Sweep.Td_Mu.mu, scalefac1*reshape(Charge_energy_trad_matrix(1,k0,:),1,ENV.Sweep.Stat.Td_Mu),'--s','linewidth',2);
    hold on;
end
xlabel('Td Mean (h)');
ylabel('Charger Energy Output (kWh)');
legend([pic,pic2],'StdEV = 10% AvgEV = 50 kWh','StdEV = 30% AvgEV = 50 kWh','StdEV = 50% AvgEV = 50 kWh',...
'StdEV = 10% AvgEV = 33 kWh','StdEV = 30% AvgEV = 33 kWh','StdEV = 50% AvgEV = 33 kWh'...
);
title('PPP');
%ylim([0.6,1]);
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

function charger_energy_mean = ExpectedChargerEnergy(ExpTrajRes, ENV, kn1, k0, k1)
    charger_energy_vec = [];
    for k2 = 1:ENV.Var_Conv.MC_trial
        for k3 = 1:ENV.Sweep.Stat.EV_MC_trial
            for k4 = 1:ENV.Sweep.Stat.Td_MC_trial
                charger_energy_vec = [charger_energy_vec, sum(ExpTrajRes{kn1,k0,k1,k2}{k3,k4}.Evdemand.value(2:end))];
            end
        end
    end
    charger_energy_mean = mean(charger_energy_vec);
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

