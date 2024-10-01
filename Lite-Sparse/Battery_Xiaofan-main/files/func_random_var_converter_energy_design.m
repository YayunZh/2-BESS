function ExpRes = func_random_var_converter_energy_design(OptAvg,CL,ENV)
%% design variance converters
% ExpRes = struct('uratio',[],'max_power',[],'process_power',[],'Conv_power_rating',[],...
%     'Bat_power_rating',[],'uratio_conv',[],'bat_eff',[],'bat_effproxy',[],'bat_eff_error',[]);

% Var_conv_energy_rating = [];
% Avg_conv_energy_rating = [];
% Bat_energy_rating = [];
% urate = [];
% max_energy = [];
% process_energy = [];
% urate_conv = [];
% bat_energy_eff = [];
% bat_energy_eff_proxy = [];
% bat_energy_eff_error = [];
% var_energy_flow = [];
% avg_energy_flow = [];
% total_energy_flow = []; 
% var_energy_process =[];
% total_energy_process = [];
% avg_energy_process = [];
% bat_output_energy_individual = [];

for k0 = 1:ENV.Sweep.Stat.Bat                     % k_0 loop for 4 battery variance
    [CL, Bat_Info_MC] = Update_k0_loop(k0,CL,ENV);       
    for k1 = 1:ENV.Sweep.Stat.Conv                % k_1 loop for 20 power converter rating
        ENV = Update_k1_loop(k1, ENV);
        for k2 = 1:ENV.Var_Conv.MC_trial       % k_2 loop is for MC simu, please
            % run('func_dc_powerflow_var_layer');
            % OutVar = func_dc_energyflow_var_layer(OptAvg,CL,ENV);
            OutVar = func_random_dc_energyflow_two_layers(ENV.Avg_Conv.e_lim_mat, Bat_Info_MC(:,k2), CL, ENV);
            
            var_energy_flow(k0,k1,k2,:,:) = OutVar.Var_energy_flow;
            avg_energy_flow(k0,k1,k2,:,:) = OutVar.Avg_energy_flow;
  %         total_energy_flow(k0,k1,k2,:,:) = OutVar.Total_energy_flow;
            
            var_energy_process(k0,k1,k2) = OutVar.Var_energy_process;
            total_energy_process(k0,k1,k2) = OutVar.Total_energy_process;
            avg_energy_process(k0,k1,k2) = OutVar.Avg_energy_process;
            
            bat_output_energy_individual(k0,k1,k2,:) = OutVar.BatOutputEnergy;
            var_conv_energy_rating(k0,k1,k2,:) = OutVar.Var_conv_energy_rating;
            avg_conv_energy_rating(k0,k1,k2,:) = ENV.Avg_Conv.e_lim_vec;
            bat_energy_rating(k0,k1,k2,:) = OutVar.Bat_energy_rating;
%             bat_energy_eff(k0,k1,k2) = OutVar.EnergyEff;
%             bat_energy_effproxy(k0,k1,k2) = OutVar.EnergyEffProxy;
%             bat_energy_eff_error(k0,k1,k2) = abs(OutVar.EnergyEff-OutVar.EnergyEffProxy);
            max_output_energy(k0,k1,k2) = abs(OutVar.Maximum_output_energy);
            bat_uratio_energy(k0,k1,k2) = abs(OutVar.Maximum_output_energy)/sum(OutVar.Bat_energy_rating);
           % process_energy(k0,k1,k2) = OutVar.Var_energy_process;
            var_conv_uratio_energy(k0,k1,k2) = OutVar.Var_energy_process/sum(OutVar.Var_conv_energy_rating);
        end 
    end
end

ExpRes.var_energy_flow = var_energy_flow;
ExpRes.avg_energy_flow = avg_energy_flow ;
%ExpRes.total_energy_flow = total_energy_flow;
ExpRes.var_energy_process = var_energy_process;
ExpRes.total_energy_process = total_energy_process;
ExpRes.avg_energy_process = avg_energy_process;
ExpRes.bat_output_energy_individual = bat_output_energy_individual;
ExpRes.var_conv_energy_rating = var_conv_energy_rating;
ExpRes.avg_conv_energy_rating = avg_conv_energy_rating;
ExpRes.bat_energy_rating = bat_energy_rating;
ExpRes.bat_uratio_energy = bat_uratio_energy;
ExpRes.max_output_energy = max_output_energy;
%ExpRes.process_energy = process_energy;
ExpRes.var_conv_uratio_energy = var_conv_uratio_energy;
% ExpRes.bat_energy_eff = bat_energy_eff;
% ExpRes.bat_energy_effproxy = bat_energy_effproxy;
% ExpRes.bat_energy_eff_error = bat_energy_eff_error;

function [CL, Bat_Info_MC] = Update_k0_loop(k0, CL, ENV)
    for i = 1:CL.Stat.Bat_num
        CL.Bat{i}.qlim_var = ENV.Sweep.Bat{i}.qlim_var(k0);
        CL.Bat{i}.qlim_mu = ENV.Sweep.Bat{i}.qlim_mu(k0);
        pd = makedist('Normal','mu',ENV.Sweep.Bat{i}.qlim_mu(k0),'sigma',ENV.Sweep.Bat{i}.qlim_var(k0));
        t = truncate(pd,0,inf);
        Bat_Info_MC(i,:) = random(t,1,ENV.Var_Conv.MC_trial);
%         Bat_Info_MC(i,:) = normrnd(CL.Bat{i}.qlim_mu, ...
%              CL.Bat{i}.qlim_var,1,ENV.Var_Conv.MC_trial);
    end 
    for j = 1:ENV.Var_Conv.MC_trial
      %  Bat_Info_MC(:,j) = sort(Bat_Info_MC(:,j));
    end
end

function ENV = Update_k1_loop(k1, ENV)
    ENV.Var_Conv.e_lim = ENV.Sweep.Conv.e_lim(k1);
end


end