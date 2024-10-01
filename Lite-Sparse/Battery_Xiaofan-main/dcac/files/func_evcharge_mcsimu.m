function [ExpTrajRes, ExpTrajRes_trad] = func_evcharge_mcsimu(CL, ENV, ENV_trad)
%% Parameters
[CL, Bat_Info_MC] = Prepare_k2_loop(CL,ENV);           % prepare the heterogenous battery samples for k2-loop to sweep
%ExpTrajRes = cell(ENV.Sweep.Stat.EV_Var, ENV.Sweep.Stat.Td_Mu, ENV.Sweep.Stat.EV_MC_trial, ENV.Sweep.Stat.Td_MC_trial);
for kn1 = 1:ENV.Sweep.Stat.EV_Mu  
    for k0 = 1:ENV.Sweep.Stat.EV_Var                   % k_0 loop for 2 ev-charging variance 
        EV_demand_Info = Update_k0_loop(kn1,k0,ENV);       
        for k1 = 1:ENV.Sweep.Stat.Td_Mu                % k_1 loop for 2 ev-charging delay
            T_delay_Info = Update_k1_loop(k1,ENV);
            parfor k2 = 1:ENV.Var_Conv.MC_trial          % k_2 loop is for MC simu of heterogenous batteries
                OutVar = func_dc_energyflow_two_layers(ENV.Avg_Conv.e_lim_mat, Bat_Info_MC(:,k2), CL, ENV); 
                OutVar_trad = func_dc_energyflow_two_layers(ENV_trad.Avg_Conv.e_lim_mat, Bat_Info_MC(:,k2), CL, ENV);
                [ExpTrajRes{kn1,k0, k1, k2}, ExpTrajRes_trad{kn1,k0, k1, k2}] = func_evcharge_mcsimu_sub(ENV, OutVar, ENV_trad, OutVar_trad, EV_demand_Info, T_delay_Info);
                % note that we use ENV instead of ENV_trad because the charging scheme are the same           
            end 
        end
    end
end

function [CL, Bat_Info_MC] = Prepare_k2_loop(CL, ENV)
    for i = 1:CL.Stat.Bat_num
        CL.Bat{i}.qlim_var = ENV.Bat.qlim_var;
        CL.Bat{i}.qlim_mu = ENV.Bat.qlim_mu;
        pd = makedist('Normal','mu',ENV.Bat.qlim_mu,'sigma',ENV.Bat.qlim_var);
        t = truncate(pd,0,inf);
        Bat_Info_MC(i,:) = random(t,1,ENV.Var_Conv.MC_trial);
%         Bat_Info_MC(i,:) = normrnd(CL.Bat{i}.qlim_mu, ...
%              CL.Bat{i}.qlim_var,1,ENV.Var_Conv.MC_trial);
    end 
    for j = 1:ENV.Var_Conv.MC_trial
        Bat_Info_MC(:,j) = sort(Bat_Info_MC(:,j));
    end
end

function EV_demand_Info = Update_k0_loop(kn1, k0, ENV)
    if (ENV.Stat.FixedDemandTdelay == 1)
        EV_demand_Info = ENV.EV_Demand;
    else
        pd = makedist('Normal','mu',ENV.Sweep.EV_Var.e_mu(kn1),'sigma',ENV.Sweep.EV_Var.e_mu(kn1)*ENV.Sweep.EV_Var.e_var_normalized(k0));
        t = truncate(pd,0.1*ENV.Sweep.EV_Var.e_mu(kn1),3*ENV.Sweep.EV_Var.e_mu(kn1));  
        EV_demand_Info = random(t,ENV.Sweep.Stat.EV_MC_trial, ENV.Stat.EV_demand_length);
    end
end

function T_delay_Info = Update_k1_loop(k1, ENV)
    if (ENV.Stat.FixedDemandTdelay == 1)
        T_delay_Info = ENV.T_delay;
    else
        pd = makedist('Exponential','mu',ENV.Sweep.Td_Mu.mu(k1));
        t = truncate(pd,0.1*ENV.Sweep.Td_Mu.mu(k1),6*ENV.Sweep.Td_Mu.mu(k1));  
        T_delay_Info = random(t,ENV.Sweep.Stat.Td_MC_trial, ENV.Stat.T_delay_length);
    end
end



end