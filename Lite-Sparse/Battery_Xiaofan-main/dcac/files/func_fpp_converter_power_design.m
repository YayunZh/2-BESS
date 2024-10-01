function ExpRes = func_fpp_converter_power_design(CL,ENV)
%% design variance converters
% ExpRes = struct('uratio',[],'max_power',[],'process_power',[],'Conv_power_rating',[],...
%     'Bat_power_rating',[],'uratio_conv',[],'bat_eff',[],'bat_effproxy',[],'bat_eff_error',[]);
% Conv_power_rating = [];
% Bat_power_rating = [];
% uratio = [];
% max_power = [];
% process_power = [];
% urate_conv = [];
% bat_eff = [];
% bat_effproxy = [];
% bat_eff_error = [];

for k0 = 1:ENV.Sweep.Stat.Bat                  % k_0 loop for 4 battery variance
    [CL, Bat_Info_MC] = Update_k0_loop(k0,CL,ENV);       
    for k1 = 1:ENV.Sweep.Stat.Conv             % k_1 loop for 20 power converter rating
        ENV = Update_k1_loop(k1, ENV);
        parfor k2 = 1:ENV.Var_Conv.MC_trial       % k_2 loop is for MC simu, please
            % run('func_dc_powerflow_var_layer');
            % OutVar = func_dc_powerflow_var_layer(OptAvg,CL,ENV);
            OutVar = func_dc_powerflow_fpp(ENV.Avg_Conv.p_lim_mat, Bat_Info_MC(:,k2), CL, ENV);
            
            var_power_flow(k0,k1,k2,:,:) = OutVar.Var_power_flow;
            avg_power_flow(k0,k1,k2,:,:) = OutVar.Avg_power_flow;
  %         total_power_flow(k0,k1,k2,:,:) = OutVar.Total_power_flow;
            
            var_power_process(k0,k1,k2) = OutVar.Var_power_process;
            total_power_process(k0,k1,k2) = OutVar.Total_power_process;
            avg_power_process(k0,k1,k2) = OutVar.Avg_power_process;
            
            bat_output_power_individual(k0,k1,k2,:) = OutVar.BatOutputpower;
            var_conv_power_rating(k0,k1,k2,:) = OutVar.Var_conv_power_rating;
            avg_conv_power_rating(k0,k1,k2,:) = ENV.Avg_Conv.p_lim_vec;
            bat_power_rating(k0,k1,k2,:) = OutVar.Bat_power_rating;
%             bat_power_eff(k0,k1,k2) = OutVar.powerEff;
%             bat_power_effproxy(k0,k1,k2) = OutVar.powerEffProxy;
%             bat_power_eff_error(k0,k1,k2) = abs(OutVar.powerEff-OutVar.powerEffProxy);
            max_output_power(k0,k1,k2) = abs(OutVar.Maximum_output_power);
            bat_uratio_power(k0,k1,k2) = abs(OutVar.Maximum_output_power)/sum(OutVar.Bat_power_rating);
           % process_power(k0,k1,k2) = OutVar.Var_power_process;
            var_conv_uratio_power(k0,k1,k2) = OutVar.Var_power_process/sum(OutVar.Var_conv_power_rating);
        end 
    end
end

ExpRes.var_power_flow = var_power_flow;
ExpRes.avg_power_flow = avg_power_flow ;
%ExpRes.total_power_flow = total_power_flow;
ExpRes.var_power_process = var_power_process;
ExpRes.total_power_process = total_power_process;
ExpRes.avg_power_process = avg_power_process;
ExpRes.bat_output_power_individual = bat_output_power_individual;
ExpRes.var_conv_power_rating = var_conv_power_rating;
ExpRes.avg_conv_power_rating = avg_conv_power_rating;
ExpRes.bat_power_rating = bat_power_rating;
ExpRes.bat_uratio_power = bat_uratio_power;
ExpRes.max_output_power = max_output_power;
%ExpRes.process_power = process_power;
ExpRes.var_conv_uratio_power = var_conv_uratio_power;
% ExpRes.bat_power_eff = bat_power_eff;
% ExpRes.bat_power_effproxy = bat_power_effproxy;
% ExpRes.bat_power_eff_error = bat_power_eff_error;


function [CL, Bat_Info_MC] = Update_k0_loop(k0, CL, ENV)
    for i = 1:CL.Stat.Bat_num
        CL.Bat{i}.curlim_var = ENV.Sweep.Bat{i}.curlim_var(k0);
        CL.Bat{i}.curlim_mu = ENV.Sweep.Bat{i}.curlim_mu(k0);
%         Bat_Info_MC(i,:) = normrnd(CL.Bat{i}.curlim_mu, ...
%              CL.Bat{i}.curlim_var,1,ENV.Var_Conv.MC_trial);
        pd = makedist('Normal','mu',ENV.Sweep.Bat{i}.curlim_mu(k0),'sigma',ENV.Sweep.Bat{i}.curlim_var(k0));
        t = truncate(pd,0,inf);
        Bat_Info_MC(i,:) = random(t,1,ENV.Var_Conv.MC_trial); 
    end 
    for j = 1:ENV.Var_Conv.MC_trial
        Bat_Info_MC(:,j) = sort(Bat_Info_MC(:,j));
    end
end

function ENV = Update_k1_loop(k1, ENV)
    ENV.Fpp_Conv.p_lim_singlevar = ENV.Sweep.Conv.p_lim(k1);
end

% function y = Update_k0_loop(k0, CL, ENV)
%     for i = 1:CL.Stat.Bat_num
%         CL.Bat{i}.curlim_var = ENV.Sweep.Bat{i}.curlim_var(k0);
%         CL.Bat{i}.curlim_mu = ENV.Sweep.Bat{i}.curlim_mu(k0);
%         CL.Bat{i}.res_var = ENV.Sweep.Bat{i}.res_var(k0);
%         CL.Bat{i}.res_mu = ENV.Sweep.Bat{i}.res_mu(k0);
%     end
%     y = CL;
% end
% 
% function y = Update_k1_loop(k1,CL,ENV) 
%     ENV.Var_Conv.p_lim = ENV.Sweep.Conv.p_lim(k1);
%     y = ENV;
% end
end