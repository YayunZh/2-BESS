function TstBat = func_battery_testing_power_v2(CL,ENV)
%% design variance converters
% TstBat = struct('urate',[],'max_power',[],'process_power',[],'Conv_power_rate',[],...
%     'Bat_power_rate',[],'urate_conv',[]);
% Conv_power_rate = [];
% Bat_power_rate = [];
% urate = [];
% max_power = [];
% process_power = [];
% urate_conv = [];
for k0 = 1:ENV.Sweep.Stat.BatPopVar
    Bat_pop_Info_MC = Update_k0_loop_2(k0, CL, ENV);
    for k05 = 1:ENV.BatPopVar.MC_trial
        CL = Update_k05_loop_2(k05, Bat_pop_Info_MC, CL, ENV);
        IniDes = func_dc_powerflow_two_layers(ENV.Avg_Conv.p_lim_mat,Bat_pop_Info_MC(:,k05),CL,ENV);
        fprintf('Given inf testing efforts, the battery power utilization is\n');
        TstBat.full_info{k0,k05}.bat_uratio = -IniDes.Maximum_output_power/sum(IniDes.Bat_power_rating);
        TstBat.full_info{k0,k05}.bat_uratio
%         ENV.Avg_Conv.p_lim_mat = IniDes.P_diff_rating_avg;
        TstBat.full_info{k0,k05}.DesignRes = IniDes;    
        TstBat.full_info{k0,k05}.SampledCL = CL;  
    %     IniDes = func_average_converter_power_design(CL,ENV);
    %     fprintf('Given inf testing efforts, with the avg layer only, the battery utilization is\n');
    %     TstBat.full_info{k0}.Avg_only.bat_uratio =  IniDes.Bat_uratio;
    %     TstBat.full_info{k0}.Avg_only.bat_uratio
    %     fprintf('Given inf testing efforts, with the avg layer only, the avg converter utilization is\n');
    %     TstBat.full_info{k0}.Avg_only.Avg_conv_uratio = IniDes.Avg_conv_uratio ;
    %     TstBat.full_info{k0}.Avg_only.Avg_conv_uratio
    %     ENV.Avg_Conv.p_lim_mat = IniDes.Conv_power_rating_partition_mat;
    %     TstBat.full_info{k0}.Avg_only.DesignRes = IniDes;
    %     
    %     SecDes = func_dc_powerflow_two_layers(CL,ENV);
    %     fprintf('Given inf testing efforts, with the two layers, the battery utilization is\n');
    %     TstBat.full_info{k0}.Two_layers.bat_uratio =  SecDes.Maximum_output_power/sum(SecDes.Bat_power_rating);
    %     TstBat.full_info{k0}.Two_layers.bat_uratio
    %     TstBat.full_info{k0}.Two_layers.DesignRes = SecDes;
        for k1 = 1:ENV.Sweep.Stat.TstBatT                  % k_1 loop for multiple battery testing cost sigma_unc    
            [CL,Bat_Info_MC] = Update_k1_loop_2(k1, CL, ENV);
            parfor k2 = 1:ENV.Tst_Bat.MC_trial            % k_2 loop for multiple battery stat MC simu
                % Redesign the Avg and Var Converters based on new battery stat
                % OutVar = func_dc_powerflow_avg_layer(CL,ENV);
                OutVar = func_dc_powerflow_two_layers(ENV.Avg_Conv.p_lim_mat,Bat_Info_MC(:,k2),CL,ENV);

                var_power_flow(k0,k05,k1,k2,:,:) = OutVar.Var_power_flow;
                avg_power_flow(k0,k05,k1,k2,:,:) = OutVar.Avg_power_flow;
      %          total_power_flow(k0,k1,k2,:,:) = OutVar.Total_power_flow;

                var_power_process(k0,k05,k1,k2) = OutVar.Var_power_process;
                total_power_process(k0,k05,k1,k2) = OutVar.Total_power_process;
                avg_power_process(k0,k05,k1,k2) = OutVar.Avg_power_process;

                bat_output_power_individual(k0,k05,k1,k2,:) = OutVar.BatOutputpower;
                var_conv_power_rating(k0,k05,k1,k2,:) = OutVar.Var_conv_power_rating;
                avg_conv_power_rating(k0,k05,k1,k2,:) = ENV.Avg_Conv.p_lim;
                bat_power_rating(k0,k05,k1,k2,:) = OutVar.Bat_power_rating;
    %             bat_power_eff(k0,k1,k2) = OutVar.powerEff;
    %             bat_power_effproxy(k0,k1,k2) = OutVar.powerEffProxy;
    %             bat_power_eff_error(k0,k1,k2) = abs(OutVar.powerEff-OutVar.powerEffProxy);
                max_output_power(k0,k05,k1,k2) = abs(OutVar.Maximum_output_power);
                bat_uratio_power(k0,k05,k1,k2) = abs(OutVar.Maximum_output_power)/sum(OutVar.Bat_power_rating);
               % process_power(k0,k1,k2) = OutVar.Var_power_process;
                var_conv_uratio_power(k0,k05,k1,k2) = OutVar.Var_power_process/sum(OutVar.Var_conv_power_rating);
            end
        end 
    end
end

TstBat.var_power_flow = var_power_flow;
TstBat.avg_power_flow = avg_power_flow ;
%ExpRes.total_power_flow = total_power_flow;
TstBat.var_power_process = var_power_process;
TstBat.total_power_process = total_power_process;
TstBat.avg_power_process = avg_power_process;
TstBat.bat_output_power_individual = bat_output_power_individual;
TstBat.var_conv_power_rating = var_conv_power_rating;
TstBat.avg_conv_power_rating = avg_conv_power_rating;
TstBat.bat_power_rating = bat_power_rating;
TstBat.bat_uratio_power = bat_uratio_power;
TstBat.max_output_power = max_output_power;
%ExpRes.process_power = process_power;
TstBat.var_conv_uratio_power = var_conv_uratio_power;
% ExpRes.bat_power_eff = bat_power_eff;
% ExpRes.bat_power_effproxy = bat_power_effproxy;
% ExpRes.bat_power_eff_error = bat_power_eff_error;


% TstBat.Conv_power_rate = Conv_power_rate;
% TstBat.Bat_power_rate = Bat_power_rate;
% TstBat.urate = urate;
% TstBat.max_power = max_power;
% TstBat.process_power = process_power;
% TstBat.urate_conv = urate_conv;

function Bat_Info_Update = Update_k0_loop_2(k0, CL, ENV2)
   for i = 1:CL.Stat.Bat_num
         Bat_Info_Update(i,:) = normrnd(ENV2.Sweep.Bat{1}.curlim_mu(k0), ...
             ENV2.Sweep.Bat{1}.curlim_pop_var(k0),1,ENV2.BatPopVar.MC_trial);
    end 
    for i = 1:ENV2.BatPopVar.MC_trial
        Bat_Info_Update(:,i) = sort(Bat_Info_Update(:,i));
    end
end

function y = Update_k05_loop_2(k05, Bat_pop_Info_MC, CL, ENV2)
    for i = 1:CL.Stat.Bat_num
%         temp = ENV2.Sweep.Bat{i}.curlim_pop_var(k0)*randn(1,1) + ENV2.Sweep.Bat{i}.curlim_mu(k0);
%         while(temp <= ENV2.Sweep.Bat{i}.curlim_mu(k0) - 3*ENV2.Sweep.Bat{i}.curlim_pop_var(k0)*randn(1,1))  % throw the data outside of 3 sigma away
%             temp = ENV2.Sweep.Bat{i}.curlim_pop_var(k0)*randn(1,1) + ENV2.Sweep.Bat{i}.curlim_mu(k0);
%         end
        CL.Bat{i}.curlim_mu = Bat_pop_Info_MC(i,k05);
        CL.Bat{i}.curlim = CL.Bat{i}.curlim_mu;
    end
    y = CL;
end


function [CL,Bat_Info_MC] = Update_k1_loop_2(k1, CL, ENV2)
    for i = 1:CL.Stat.Bat_num
        CL.Bat{i}.curlim_var = ENV2.Sweep.Bat{i}.curlim_test_var(k1);
        %Bat_stat(i,:) = normrnd(CL.Bat{i}.curlim_mu, ENV2.Sweep.TstBatT(k1)*CL.Bat{i}.curlim_var, 1, ENV2.Tst_Bat.MC_trial,1);
    end  % we generate the battery statas together, so that each converter design handles the same battery sets
    for i = 1:CL.Stat.Bat_num
         Bat_Info_MC(i,:) = normrnd(CL.Bat{i}.curlim_mu, ...
             CL.Bat{i}.curlim_var,1,ENV2.Tst_Bat.MC_trial);
    end 
    for i = 1:ENV2.Tst_Bat.MC_trial
        Bat_Info_MC(:,i) = sort(Bat_Info_MC(:,i));
    end
end

% function y = Update_k2_loop_2(k2,CL)
%     for l = 1:CL.Stat.Bat_num
%         %CL.Bat{l}.curlim = mvnrnd([CL.Bat{l}.curlim_mu,CL.Bat{l}.volt_mu], diag([CL.Bat{l}.curlim_var,CL.Bat{l}.volt_var]),1);
%          CL.Bat{l}.curlim = normrnd(CL.Bat{l}.curlim_mu, CL.Bat{l}.curlim_var);
%     end
%     CL.Cap{1}.volt = CL.Bus{1}.volt - CL.Bat{1}.volt - CL.Bat{2}.volt - CL.Bat{3}.volt...
%     - CL.Bat{4}.volt - CL.Bat{5}.volt - CL.Bat{6}.volt...
%     - CL.Bat{7}.volt - CL.Bat{8}.volt - CL.Bat{9}.volt;     % Decide the capactior votlages
%     y = CL;
% end

% function y = merge_multilayer_flow(Avgc,Varc,CL)
%     Avgc(CL.Bus{1}.ind:end,:) = 0;
%     Avgc(:,CL.Bus{1}.ind:end) = 0;
%     y = Avgc + Varc;
% end

end