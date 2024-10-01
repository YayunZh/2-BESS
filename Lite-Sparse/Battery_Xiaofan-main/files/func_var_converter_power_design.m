function ExpRes = func_var_converter_power_design(CL,ENV)
%% design variance converters
% ExpRes = struct('uratio',[],'max_power',[],'process_power',[],'Conv_power_rating',[],...
%     'Bat_power_rating',[],'uratio_conv',[],'bat_eff',[],'bat_effproxy',[],'bat_eff_error',[]);

k0_his = 1;   % k0_his selects different battery power variance, k0_his = 1 reprents the least variance
k2_his = 1;   % k2_his loop selects different battery stat and power converter rating stat

conv_power_rating = [];
bat_power_rating = [];
urate = [];
max_power = [];
process_power = [];
urate_conv = [];
bat_power_eff = [];
bat_power_effproxy = [];
bat_power_eff_error = [];

for k0 = 1:ENV.Sweep.Stat.Loss             % k_0 loop for different efficiency variance
    CL = Update_k0_loop(k0,CL,ENV);       
    for k1 = 1:ENV.Sweep.Stat.Conv        % k_1 loop for different power converter rating
        for k2 = 1:ENV.PostEff.MC_trial   % k_2 loop is for MC simu, please
            %run('func_dc_powerflow_var_layer')
            %ExpRes.Var_power_flow(k0,k1,k2) = OutVar.Var_power_flow;
            CL = Update_k2_loop(k2,CL,ENV);
%            conv_power_rating(k0,k1,k2,:) = ExpRes.conv_power_rating(k0_his,k1,k2_his,:);
%            bat_power_rating(k0,k1,k2,:) = ExpRes.bat_power_rating(k0_his,k1,k2_his,:);
            for j = 1: CL.Stat.Bat_num
                BatpowerLoss(j) = ExpRes.bat_output_power_individual(k0_his,k1,k2_his,j)*CL.Bat{j}.ploss_coeff;
                BatpowerLossProxy(j) = ExpRes.bat_output_power_individual(k0_his,k1,k2_his,j)*CL.Bat{j}.ploss_coeff_mu;
            end
            bat_power_eff(k0,k1,k2) = 1 - sum(BatpowerLoss)/(ExpRes.max_output_power(k0_his,k1,k2_his));
            bat_power_effproxy(k0,k1,k2) = 1 - sum(BatpowerLossProxy)/(ExpRes.max_output_power(k0_his,k1,k2_his));
            bat_power_eff_error(k0,k1,k2) = abs(bat_power_eff(k0,k1,k2)-bat_power_effproxy(k0,k1,k2))/(bat_power_eff(k0,k1,k2));
            max_output_power(k0,k1,k2) = ExpRes.max_output_power(k0_his,k1,k2_his);
            bat_uratio_power(k0,k1,k2) = ExpRes.bat_uratio_power(k0_his,k1,k2_his);
        end 
    end
end

%ExpResPlot.conv_power_rating = conv_power_rating;
%ExpResPlot.bat_power_rating = bat_power_rating;
ExpResPlot.bat_power_eff = bat_power_eff;
ExpResPlot.bat_power_effproxy = bat_power_effproxy;
ExpResPlot.bat_power_eff_error = bat_power_eff_error;
ExpResPlot.max_output_power = max_power;
ExpResPlot.bat_uratio_power = bat_uratio_power;

function y = Update_k0_loop(k0, CL, ENV)
    for i = 1:CL.Stat.Bat_num
        CL.Bat{i}.ploss_coeff_var = ENV.Sweep.Bat{i}.ploss_coeff_var(k0);
        CL.Bat{i}.ploss_coeff_mu = ENV.Sweep.Bat{i}.ploss_coeff_mu(k0);
    end
    y = CL;
end

function y = Update_k2_loop(k2,CL,ENV)
    for l = 1:CL.Stat.Bat_num
        CL.Bat{l}.ploss_coeff = normrnd(CL.Bat{l}.ploss_coeff_mu, CL.Bat{l}.ploss_coeff_var,1,1);
    end
    y = CL;
end


end