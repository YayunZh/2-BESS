function DiagnosisBat = func_battery_non_diagnostic_energy(CL,ENV)
%% design variance converters
% TstBat = struct('urate',[],'max_power',[],'process_power',[],'Conv_power_rate',[],...
%     'Bat_power_rate',[],'urate_conv',[]);
% Conv_power_rate = [];
% Bat_power_rate = [];
% urate = [];
% max_power = [];
% process_power = [];
% urate_conv = [];
for k0 = 1: ENV.Sweep.Stat.BatPopVar
    for k05 = 1:ENV.BatPopVar.MC_trial
        CL = Update_k05_loop_2(k0,CL,ENV);
        IniDes = func_two_converter_energy_design(CL,ENV);
        fprintf('Given inf testing efforts, with the avg layer only, the battery utilization is\n');
        DiagnosisBat.full_info{k0,k05}.bat_uratio = IniDes.Maximum_output_energy/sum(IniDes.Bat_energy_rating);
        DiagnosisBat.full_info{k0,k05}.bat_uratio
        ENV.Avg_Conv.e_lim_mat = IniDes.P_diff_rating_avg;
        DiagnosisBat.full_info{k0,k05}.DesignRes = IniDes;    
        DiagnosisBat.full_info{k0,k05}.SampledCL = CL;  
    %     IniDes = func_average_converter_energy_design(CL,ENV);
    %     fprintf('Given inf testing efforts, with the avg layer only, the battery utilization is\n');
    %     TstBat.full_info{k0}.Avg_only.bat_uratio =  IniDes.Bat_uratio;
    %     TstBat.full_info{k0}.Avg_only.bat_uratio
    %     fprintf('Given inf testing efforts, with the avg layer only, the avg converter utilization is\n');
    %     TstBat.full_info{k0}.Avg_only.Avg_conv_uratio = IniDes.Avg_conv_uratio ;
    %     TstBat.full_info{k0}.Avg_only.Avg_conv_uratio
    %     ENV.Avg_Conv.e_lim_mat = IniDes.Conv_energy_rating_partition_mat;
    %     TstBat.full_info{k0}.Avg_only.DesignRes = IniDes;
    %     
    %     SecDes = func_dc_energyflow_two_layers(CL,ENV);
    %     fprintf('Given inf testing efforts, with the two layers, the battery utilization is\n');
    %     TstBat.full_info{k0}.Two_layers.bat_uratio =  SecDes.Maximum_output_energy/sum(SecDes.Bat_energy_rating);
    %     TstBat.full_info{k0}.Two_layers.bat_uratio
    %     TstBat.full_info{k0}.Two_layers.DesignRes = SecDes;  
        for k1= 1:ENV.Sweep.Stat.BatDiagnosisNon                   % k_1 loop for multiple battery diagnosis sigma_unc    
            [CL,Bat_Info_MC] = Update_k1_loop_2(k1, CL, ENV);
            parfor k2 = 1:ENV.Diagnosis_Bat.MC_trial            % k_2 loop for multiple battery stat MC simu
                % Redesign the Avg and Var Converters based on new battery stat
                % OutVar = func_dc_energyflow_avg_layer(CL,ENV);
                OutVar = func_dc_energyflow_check_two_layers(ENV.Avg_Conv.e_lim_mat,Bat_Info_MC(:,k2),IniDes, CL, ENV);

                var_energy_flow(k0,k05,k1,k2,:,:) = OutVar.Var_energy_flow;
                avg_energy_flow(k0,k05,k1,k2,:,:) = OutVar.Avg_energy_flow;
                %total_energy_flow(k0,k05,k1,k2,:,:) = OutVar.Total_energy_flow;
      %          total_energy_flow(k0,k1,k2,:,:) = OutVar.Total_energy_flow;

                var_energy_process(k0,k05,k1,k2) = OutVar.Var_energy_process;
                total_energy_process(k0,k05,k1,k2) = OutVar.Total_energy_process;
                avg_energy_process(k0,k05,k1,k2) = OutVar.Avg_energy_process;

                bat_output_energy_individual(k0,k05,k1,k2,:) = OutVar.BatOutputEnergy;
                var_conv_energy_rating(k0,k05,k1,k2,:) = OutVar.Var_conv_energy_rating;
                avg_conv_energy_rating(k0,k05,k1,k2,:) = ENV.Avg_Conv.e_lim;
                bat_energy_rating(k0,k05,k1,k2,:) = OutVar.Bat_energy_rating;
                
                lambda_shrink(k0,k05,k1,k2,:) = OutVar.lambda_shrink;
    %             bat_energy_eff(k0,k1,k2) = OutVar.EnergyEff;
    %             bat_energy_effproxy(k0,k1,k2) = OutVar.EnergyEffProxy;
    %             bat_energy_eff_error(k0,k1,k2) = abs(OutVar.EnergyEff-OutVar.EnergyEffProxy);
                max_output_energy(k0,k05,k1,k2) = abs(OutVar.Maximum_output_energy);
                bat_uratio_energy(k0,k05,k1,k2) = abs(OutVar.Maximum_output_energy)/sum(OutVar.Bat_energy_rating);
                if(bat_uratio_energy(k0,k05,k1,k2) >= 1)
                    fprintf('wrong battery uratio');
                    pause;
                end
               % process_energy(k0,k1,k2) = OutVar.Var_energy_process;
                var_conv_uratio_energy(k0,k05,k1,k2) = OutVar.Var_energy_process/sum(OutVar.Var_conv_energy_rating);
            end
        end 
    end
end

DiagnosisBat.var_energy_flow = var_energy_flow;
DiagnosisBat.avg_energy_flow = avg_energy_flow ;
DiagnosisBat.lambda_shrink = lambda_shrink;
%ExpRes.total_energy_flow = total_energy_flow;
DiagnosisBat.var_energy_process = var_energy_process;
DiagnosisBat.total_energy_process = total_energy_process;
DiagnosisBat.avg_energy_process = avg_energy_process;
DiagnosisBat.bat_output_energy_individual = bat_output_energy_individual;
DiagnosisBat.var_conv_energy_rating = var_conv_energy_rating;
DiagnosisBat.avg_conv_energy_rating = avg_conv_energy_rating;
DiagnosisBat.bat_energy_rating = bat_energy_rating;
DiagnosisBat.bat_uratio_energy = bat_uratio_energy;
DiagnosisBat.max_output_energy = max_output_energy;
%ExpRes.process_energy = process_energy;
DiagnosisBat.var_conv_uratio_energy = var_conv_uratio_energy;
% ExpRes.bat_energy_eff = bat_energy_eff;
% ExpRes.bat_energy_effproxy = bat_energy_effproxy;
% ExpRes.bat_energy_eff_error = bat_energy_eff_error;


% TstBat.Conv_power_rate = Conv_power_rate;
% TstBat.Bat_power_rate = Bat_power_rate;
% TstBat.urate = urate;
% TstBat.max_power = max_power;
% TstBat.process_power = process_power;
% TstBat.urate_conv = urate_conv;

function y = Update_k05_loop_2(k0, CL, ENV2)
    for i = 1:CL.Stat.Bat_num
%         temp = ENV2.Sweep.Bat{i}.qlim_pop_var(k0)*randn(1,1) + ENV2.Sweep.Bat{i}.qlim_mu(k0);
%         while(temp <= ENV2.Sweep.Bat{i}.qlim_mu(k0) - 3*ENV2.Sweep.Bat{i}.qlim_pop_var(k0)*randn(1,1))  % throw the data outside of 3 sigma away
%             temp = ENV2.Sweep.Bat{i}.qlim_pop_var(k0)*randn(1,1) + ENV2.Sweep.Bat{i}.qlim_mu(k0);
%         end
        CL.Bat{i}.qlim_mu = ENV2.Sweep.Bat{i}.qlim_pop_var(k0)*randn(1,1) + ENV2.Sweep.Bat{i}.qlim_mu(k0);
        CL.Bat{i}.qlim = CL.Bat{i}.qlim_mu;
    end
    y = CL;
end

function [CL,Bat_Info_MC] = Update_k1_loop_2(k1, CL, ENV2)
    for i = 1:CL.Stat.Bat_num
        CL.Bat{i}.qlim_var = ENV2.Sweep.Bat{i}.qlim_diagnosis_var(k1);
        %Bat_stat(i,:) = normrnd(CL.Bat{i}.curlim_mu, ENV2.Sweep.TstBatT(k1)*CL.Bat{i}.curlim_var, 1, ENV2.Tst_Bat.MC_trial,1);
    end  % we generate the battery statas together, so that each converter design handles the same battery sets
    for i = 1:CL.Stat.Bat_num
         Bat_Info_MC(i,:) = normrnd(CL.Bat{i}.qlim_mu*ENV2.Sweep.DiscountFactor(k1), CL.Bat{i}.qlim_var,1,ENV2.Diagnosis_Bat.MC_trial);
    end 
end

function y = Update_k2_loop_2(k2,CL)
    for l = 1:CL.Stat.Bat_num
        %CL.Bat{l}.qlim = mvnrnd([CL.Bat{l}.qlim_mu,CL.Bat{l}.volt_mu], diag([CL.Bat{l}.qlim_var,CL.Bat{l}.volt_var]),1);
         CL.Bat{l}.qlim = normrnd(CL.Bat{l}.qlim_mu, CL.Bat{l}.qlim_var);
    end
    CL.Cap{1}.volt = CL.Bus{1}.volt - CL.Bat{1}.volt - CL.Bat{2}.volt - CL.Bat{3}.volt...
    - CL.Bat{4}.volt - CL.Bat{5}.volt - CL.Bat{6}.volt...
    - CL.Bat{7}.volt - CL.Bat{8}.volt - CL.Bat{9}.volt;     % Decide the capactior votlages
    y = CL;
end

% function y = merge_multilayer_flow(Avgc,Varc,CL)
%     Avgc(CL.Bus{1}.ind:end,:) = 0;
%     Avgc(:,CL.Bus{1}.ind:end) = 0;
%     y = Avgc + Varc;
% end

end