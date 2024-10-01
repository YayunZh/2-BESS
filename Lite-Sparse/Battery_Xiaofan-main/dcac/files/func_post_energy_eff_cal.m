function ExpResPlot = func_post_energy_eff_cal(ExpRes,CL,ENV)
%% design variance converters
% ExpRes = struct('uratio',[],'max_power',[],'process_power',[],'Conv_power_rating',[],...
%     'Bat_power_rating',[],'uratio_conv',[],'bat_eff',[],'bat_effproxy',[],'bat_eff_error',[]);

k0_his = 1;   % k0_his = 1 selects the least variance battery
k2_his = 1;   % k2_his = 1 selects the first of the MC simulation 

conv_energy_rating = [];
bat_energy_rating = [];
urate = [];
max_energy = [];
process_energy = [];
urate_conv = [];
bat_energy_eff = [];
bat_energy_effproxy = [];
bat_energy_eff_error = [];

for k0 = 1:ENV.Sweep.Stat.Bat             % k_0 loop for different efficiency variance
    CL = Update_k0_loop(k0,CL,ENV);       
    for k1 = 1:ENV.Sweep.Stat.Conv        % k_1 loop for different power converter rating
        for k2 = 1:ENV.PostEff.MC_trial   % k_2 loop is for MC simu, please
            %run('func_dc_powerflow_var_layer')
            %ExpRes.Var_power_flow(k0,k1,k2) = OutVar.Var_power_flow;
            CL = Update_k2_loop(k2,CL,ENV);
%           conv_energy_rating(k0,k1,k2,:) = ExpRes.conv_energy_rating(k0_his,k1,k2_his,:);
%           bat_energy_rating(k0,k1,k2,:) = ExpRes.bat_energy_rating(k0_his,k1,k2_his,:);
            for j = 1: CL.Stat.Bat_num
                BatEnergyLoss(j) = ExpRes.bat_output_energy_individual(j)*CL.Bat{j}.eloss_coeff;
                BatEnergyLossProxy(j) = ExpRes.bat_output_energy_individual(j)*CL.Bat{j}.eloss_coeff_mu;
            end
            bat_energy_eff(k0,k1,k2) = 1 - sum(BatEnergyLoss)/(ExpRes.max_output_energy(k0_his,k1,k2_his));
            bat_energy_effproxy(k0,k1,k2) = 1 - sum(BatEnergyLossProxy)/(ExpRes.max_output_energy(k0_his,k1,k2_his));
            bat_energy_eff_error(k0,k1,k2) = abs(bat_energy_eff(k0,k1,k2)-bat_energy_effproxy(k0,k1,k2))/(bat_energy_eff(k0,k1,k2));
            max_output_energy(k0,k1,k2) = ExpRes.max_output_energy(k0_his,k1,k2_his);
            bat_uratio_energy(k0,k1,k2) = ExpRes.bat_uratio_energy(k0_his,k1,k2_his);
        end 
    end
end

%ExpResPlot.conv_energy_rating = conv_energy_rating;
%ExpResPlot.bat_energy_rating = bat_energy_rating;
ExpResPlot.bat_energy_eff = bat_energy_eff;
ExpResPlot.bat_energy_effproxy = bat_energy_effproxy;
ExpResPlot.bat_energy_eff_error = bat_energy_eff_error;
ExpResPlot.max_output_energy = max_energy;
ExpResPlot.bat_uratio_energy = bat_uratio_energy;

function y = Update_k0_loop(k0, CL, ENV)
    for i = 1:CL.Stat.Bat_num
        CL.Bat{i}.eloss_coeff_var = ENV.Sweep.Bat{i}.eloss_coeff_var(k0);
        CL.Bat{i}.eloss_coeff_mu = ENV.Sweep.Bat{i}.eloss_coeff_mu(k0);
    end
    y = CL;
end

function y = Update_k2_loop(k2,CL,ENV)
    for l = 1:CL.Stat.Bat_num
        CL.Bat{l}.eloss_coeff = normrnd(CL.Bat{l}.eloss_coeff_mu, CL.Bat{l}.eloss_coeff_var,1,1);
    end
    y = CL;
end
end