function y = func_plot_task4_7(ENV,ExpTrajRes,ExpTrajRes_trad)

scalefac1 = 150/36;
% ExpTrajRes_typ = ExpTrajRes{1,1,1}{1,1};
% ExpTrajRes_typ_trad = ExpTrajRes_trad{1,1,1}{1,1};

% ExpTrajRes_typ.Pbessout = ExpectedTrajPbessout (ExpTrajRes, ENV);
% ExpTrajRes_typ_trad.Pbessout = ExpectedTrajPbessout (ExpTrajRes_trad, ENV);
for kn1 = 1:ENV.Sweep.Stat.EV_Mu                       % k_n1 loop for ev demand mean      
    for k0 = 1:ENV.Sweep.Stat.EV_Var                  % k_0 loop for ev demand variance      
        for k1 = 1:ENV.Sweep.Stat.Td_Mu             % k_1 loop for Tdelay mu 
            Uenergy  = ExpectedTrajuratio_energy(ExpTrajRes, ENV, kn1, k0, k1);
            Uenergy_trad = ExpectedTrajuratio_energy(ExpTrajRes_trad, ENV, kn1, k0, k1);
            
            U_ratio_matrix(kn1,k0,k1) = Uenergy.mean;
            U_ratio_trad_matrix(kn1,k0,k1) = Uenergy_trad.mean;
            
            U_ratio_std_matrix(kn1,k0,k1) = Uenergy.std;
            U_ratio_std_trad_matrix(kn1,k0,k1) = Uenergy_trad.std;
        end
    end
end

%% Expected Battery Energy Utilization
figure();
for k1 = 1:ENV.Sweep.Stat.Td_Mu
    pic(k1) = plot(ENV.Sweep.EV_Var.e_var_normalized, U_ratio_matrix(2,:,k1),'-s','linewidth',2);
    hold on;
end
for k1 = 1:ENV.Sweep.Stat.Td_Mu
    pic2(k1) = plot(ENV.Sweep.EV_Var.e_var_normalized, U_ratio_matrix(1,:,k1),'--s','linewidth',2);
    hold on;
end
xlabel('EV Demand Deviation');
ylabel('Battery Energy Utilization');
legend([pic,pic2],'$\bar{T}_d$ = 0.5 h $\bar{E}_{ev}$ = 50 kWh',...
'$\bar{T}_d$ = 1.5 h $\bar{E}_{ev}$ = 50 kWh',... 
'$\bar{T}_d$ = 2.5 h $\bar{E}_{ev}$ = 50 kWh',...
'$\bar{T}_d$ = 0.5 h $\bar{E}_{ev}$ = 33 kWh',...
'$\bar{T}_d$ = 1.5 h $\bar{E}_{ev}$ = 33 kWh',...
'$\bar{T}_d$ = 2.5 h $\bar{E}_{ev}$ = 33 kWh','Interpreter','latex');
title('LS-HiPPP');
ylim([0.6,1]);
grid on;
grid minor;

figure();
for k1 = 1:ENV.Sweep.Stat.Td_Mu
    pic(k1) = plot(ENV.Sweep.EV_Var.e_var_normalized, U_ratio_trad_matrix(2,:,k1),'-d','linewidth',2);
    hold on;
end
for k1 = 1:ENV.Sweep.Stat.Td_Mu
    pic2(k1) = plot(ENV.Sweep.EV_Var.e_var_normalized, U_ratio_matrix(1,:,k1),'--d','linewidth',2);
    hold on;
end
xlabel('EV Demand Deviation');
ylabel('Battery Energy Utilization');
legend([pic,pic2],'$\bar{T}_d$ = 0.5 h $\bar{E}_{ev}$ = 50 kWh',...
'$\bar{T}_d$ = 1.5 h $\bar{E}_{ev}$ = 50 kWh',... 
'$\bar{T}_d$ = 2.5 h $\bar{E}_{ev}$ = 50 kWh',...
'$\bar{T}_d$ = 0.5 h $\bar{E}_{ev}$ = 33 kWh',...
'$\bar{T}_d$ = 1.5 h $\bar{E}_{ev}$ = 33 kWh',...
'$\bar{T}_d$ = 2.5 h $\bar{E}_{ev}$ = 33 kWh','Interpreter','latex');
title('PPP');
ylim([0.6,1]);
grid on;
grid minor;

%% Battery Energy Deviations
% figure();
% for k1 = 1:ENV.Sweep.Stat.Td_Mu
%     pic(k1) = plot(ENV.Sweep.EV_Var.e_var_normalized, U_ratio_std_matrix(2,:,k1),'-s','linewidth',2);
%     hold on;
% end
% for k1 = 1:ENV.Sweep.Stat.Td_Mu
%     pic2(k1) = plot(ENV.Sweep.EV_Var.e_var_normalized, U_ratio_std_matrix(1,:,k1),'--s','linewidth',2);
%     hold on;
% end
% xlabel('EV Demand Variation');
% ylabel('Expected Battery Energy Utilization');
% legend([pic,pic2],'AvgTd = 0.5h AvgEV = 50 kWh','AvgTd = 1.5h AvgEV = 50 kWh','AvgTd = 2.5h AvgEV = 50 kWh',...
% 'AvgTd = 0.5h AvgEV = 33 kWh','AvgTd = 1.5h AvgEV = 33 kWh','AvgTd = 2.5h AvgEV = 33 kWh'...
% );
% title('LS-HiPPP');
% ylim([0.6,1]);
% grid on;
% grid minor;
% 
% figure();
% for k1 = 1:ENV.Sweep.Stat.Td_Mu
%     pic(k1) = plot(ENV.Sweep.EV_Var.e_var_normalized, U_ratio_std_trad_matrix(2,:,k1),'-d','linewidth',2);
%     hold on;
% end
% for k1 = 1:ENV.Sweep.Stat.Td_Mu
%     pic2(k1) = plot(ENV.Sweep.EV_Var.e_var_normalized, U_ratio_std_matrix(1,:,k1),'--d','linewidth',2);
%     hold on;
% end
% xlabel('EV Demand Variation');
% ylabel('Expected Battery Energy Utilization');
% legend([pic,pic2],'AvgTd = 0.5h AvgEV = 50 kWh','AvgTd = 1.5h AvgEV = 50 kWh','AvgTd = 2.5h AvgEV = 50 kWh',...
% 'AvgTd = 0.5h AvgEV = 33 kWh','AvgTd = 1.5h AvgEV = 33 kWh','AvgTd = 2.5h AvgEV = 33 kWh'...
% );
% title('PPP');
% ylim([0.6,1]);
% grid on;
% grid minor;


%% User Defined Functions
function Uenergy = ExpectedTrajuratio_energy(ExpTrajRes, ENV, kn1, k0, k1)
    uratio_energy_vec = [];
    for k2 = 1:ENV.Var_Conv.MC_trial
        for k3 = 1:ENV.Sweep.Stat.EV_MC_trial
            for k4 = 1:ENV.Sweep.Stat.Td_MC_trial
                uratio_energy_vec = [uratio_energy_vec, ExpTrajRes{kn1,k0,k1,k2}{k3,k4}.uratio_energy.value(2:end)];
            end
        end
    end
    Uenergy.mean = mean(uratio_energy_vec);
    Uenergy.std = std(uratio_energy_vec);
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

