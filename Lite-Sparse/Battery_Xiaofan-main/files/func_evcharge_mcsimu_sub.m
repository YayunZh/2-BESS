function [ExpTrajRes, ExpTrajRes_trad] = func_evcharge_mcsimu_sub(ENV, OutVar, ENV_trad, OutVar_trad, EV_demand_Info, T_delay_Info)
    for k3 = 1:ENV.Sweep.Stat.EV_MC_trial
        for k4 = 1:ENV.Sweep.Stat.Td_MC_trial
            ExpTrajRes{k3, k4} = func_EVcharging(OutVar, ENV, EV_demand_Info(k3,:), T_delay_Info(k4,:));
            ExpTrajRes_trad{k3, k4} = func_EVcharging(OutVar_trad, ENV, EV_demand_Info(k3,:), T_delay_Info(k4,:));  % note that we use ENV instead of ENV_trad because the charging scheme a
    %                      ExpTrajRes{k0, k1, k2} = func_EVcharging(OutVar, ENV, EV_demand_Info(1,:), T_delay_Info(1,:));
    %                      ExpTrajRes_trad{k0, k1, k2} = func_EVcharging(OutVar_trad, ENV, EV_demand_Info(1,:), T_delay_Info(1,:));
        end
    end 
end