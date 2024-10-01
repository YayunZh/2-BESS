function ExpRes = func_lsv2_converter_energy_design(CL,ENV,bool_sort)

% battery creation
% CL = Flatten_battery_stat(CL); % use the flatened distributions of q
% CL = Update_cap_volt(CL);

for k0 = 1:ENV.Sweep.Stat.Bat                  % k0 loop for 8 battery variance, n of bat. types
    [CL, Bat_Info_MC] = Update_k0_loop(k0,CL,ENV,bool_sort);    % 9*100 for each k(bat. type), qlim mu var 
    for k1 = 1:ENV.Sweep.Stat.Conv             % k_1 loop for 20 energy converter rating
        [ENV,err] = Update_k1_loop(k1, CL,ENV,Bat_Info_MC); %e_lim update, 1 * 100
        if (err)
            continue; %if for ko bet variance, k1 converter rating limit is invalid, skip
        end
        parfor k2 = 1:ENV.Fir_Conv.MC_trial       % k_2 loop is for MC simu, please
            % run('func_dc_energyflow_var_layer');
            % OutVar = func_dc_energyflow_var_layer(OptAvg,CL,ENV);
            
           
            OutVar = func_dc_energyflow_lsv2(k2,Bat_Info_MC(:,k2), CL, ENV);
            
            % var_energy_flow(k0,k1,k2,:,:) = OutVar.Var_energy_flow;
            avg_energy_flow(k0,k1,k2,:,:) = OutVar.Avg_energy_flow;
            lsv_energy_flow(k0,k1,k2,:,:) = OutVar.Lsv_energy_flow;
  %         total_energy_flow(k0,k1,k2,:,:) = OutVar.Total_energy_flow;
            
            % var_energy_process(k0,k1,k2) = OutVar.Var_energy_process;
            total_energy_process(k0,k1,k2) = OutVar.Total_energy_process;
            lsv_energy_process(k0,k1,k2) = OutVar.Lsv_energy_process;

            avg_energy_process(k0,k1,k2) = OutVar.Avg_energy_process;
            
            bat_output_energy_individual(k0,k1,k2,:) = OutVar.BatOutputEnergy;
            % var_conv_energy_rating(k0,k1,k2,:) = OutVar.Var_conv_energy_rating;
            % avg_conv_energy_rating(k0,k1,k2,:) = ENV.Avg_Conv.e_lim_vec;
            bat_energy_rating(k0,k1,k2,:) = OutVar.Bat_energy_rating;
%             bat_energy_eff(k,k1,k2) = OutVar.energyEff;
%             bat_energy_effproxy(k,k1,k2) = OutVar.energyEffProxy;
%             bat_energy_eff_error(k,k1,k2) = abs(OutVar.energyEff-OutVar.energyEffProxy);
            max_output_energy(k0,k1,k2) = abs(OutVar.Maximum_output_energy);
            bat_uratio_energy(k0,k1,k2) = abs(OutVar.Maximum_output_energy)/sum(OutVar.Bat_energy_rating);
           % process_energy(k,k1,k2) = OutVar.Var_energy_process;
           % var_conv_uratio_energy(k0,k1,k2) = OutVar.Var_energy_process/sum(OutVar.Var_conv_energy_rating);
        end 
    end

end

% ExpRes.var_energy_flow = var_energy_flow;
% ExpRes.avg_energy_flow = avg_energy_flow ;
% %ExpRes.total_energy_flow = total_energy_flow;
% ExpRes.var_energy_process = var_energy_process;
ExpRes.total_energy_process = total_energy_process;
ExpRes.lsv_energy_process = lsv_energy_process;
ExpRes.bat_output_energy_individual = bat_output_energy_individual;
% ExpRes.var_conv_energy_rating = var_conv_energy_rating;
% ExpRes.avg_conv_energy_rating = avg_conv_energy_rating;
ExpRes.bat_energy_rating = bat_energy_rating;
ExpRes.bat_uratio_energy = bat_uratio_energy;
ExpRes.max_output_energy = max_output_energy;
%ExpRes.process_energy = process_energy;
% ExpRes.var_conv_uratio_energy = var_conv_uratio_energy;
% ExpRes.bat_energy_eff = bat_energy_eff;
% ExpRes.bat_energy_effproxy = bat_energy_effproxy;
% ExpRes.bat_energy_eff_error = bat_energy_eff_error;
end

function [CL, Bat_Info_MC] = Update_k0_loop(k0, CL, ENV,bool_sort)
    for i = 1:CL.Stat.Bat_num % 9
        CL.Bat{i}.qlim_var = ENV.Sweep.Bat{i}.qlim_var(k0);
        CL.Bat{i}.qlim_mu = ENV.Sweep.Bat{i}.qlim_mu(k0);
        
        pd = makedist('Normal','mu',ENV.Sweep.Bat{i}.qlim_mu(k0),'sigma',ENV.Sweep.Bat{i}.qlim_var(k0));
        t = truncate(pd,0,inf);
        Bat_Info_MC(i,:) = random(t,1,ENV.Fir_Conv.MC_trial);
%         Bat_Info_MC(i,:) = normrnd(CL.Bat{i}.qlim_mu, ...
%              CL.Bat{i}.qlim_var,1,ENV.Var_Conv.MC_trial);
    end 
    if (bool_sort)
        for j = 1:ENV.Fir_Conv.MC_trial
            Bat_Info_MC(:,j) = sort(Bat_Info_MC(:,j));
        end
    end
end

function [ENV,err] = Update_k1_loop(k1, CL, ENV,Bat_Info_MC)
    % max_op_energy = zeros(4,ENV.Fir_Conv.MC_trial);
    % max_op_energy(1,:) = CL.Bat{1}.volt .* Bat_Info_MC(1,:);
    % max_op_energy(2,:) = CL.Bat{2}.volt .* Bat_Info_MC(2,:);
    % max_op_energy(3,:)= CL.Bat{8}.volt .* Bat_Info_MC(8,:);
    % max_op_energy(4,:) = CL.Bat{9}.volt .* Bat_Info_MC(9,:);
    % avg_energy = 0.25 .* (sum(max_op_energy)); % sum of each column
    % ENV.Sec_Conv.e_lim_single = max(abs(max_op_energy - avg_energy));
    % %%%%%%%%%%% 
    % for i = 1:100
    %     if (ENV.Sweep.Conv_energy_sum(k1) < 4 * ENV.Sec_Conv.e_lim_single(i))
    %         ENV.Sec_Conv.e_lim_single(i)  = ENV.Sweep.Conv_energy_sum(k1)/4;
    %         ENV.Fir_Conv.e_lim_single(i) = 0;
    %     else
    %         ENV.Fir_Conv.e_lim_single(i) = (ENV.Sweep.Conv_energy_sum(k1) - 4 * ENV.Sec_Conv.e_lim_single(i))/ENV.Fir_Conv.Num; % 8
    %     end
    % end

    max_op_energy = zeros(4,1);
    max_op_energy(1) = CL.Bat{1}.volt .* mean(Bat_Info_MC(1,:));
    max_op_energy(2) = CL.Bat{2}.volt .* mean(Bat_Info_MC(2,:));
    max_op_energy(3)= CL.Bat{8}.volt .* mean(Bat_Info_MC(8,:));
    max_op_energy(4) = CL.Bat{9}.volt .* mean(Bat_Info_MC(9,:));
    avg_energy = 0.25 .* (sum(max_op_energy)); % sum of each column
    ENV.Sec_Conv.e_lim_single = zeros(1,ENV.Fir_Conv.MC_trial);
    ENV.Sec_Conv.e_lim_single(:) = 1 * max(abs(max_op_energy - avg_energy));
    %%%%%%%%%%% 
    ENV.Fir_Conv.e_lim_single = zeros(1,ENV.Fir_Conv.MC_trial);
    ENV.Fir_Conv.e_lim_single(:) = (ENV.Sweep.Conv_energy_sum(k1) - 4 * ENV.Sec_Conv.e_lim_single(1))/ENV.Fir_Conv.Num; % 
    err = 0;
    if (ENV.Fir_Conv.e_lim_single(1) < 0)
            err = 1;
            ENV.Sec_Conv.e_lim_single(:) = 0;
            ENV.Fir_Conv.e_lim_single(:) = ENV.Sweep.Conv_energy_sum(k1)./ENV.Fir_Conv.Num;
     end
end