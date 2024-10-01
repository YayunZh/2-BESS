function OutRes = func_two_converter_energy_design(CL,ENV)
%Input the diffferetial energy delivery connection
trial_num = ENV.Avg_Conv.trial_num;
converter_num = ENV.Avg_Conv.Num;

comb_explosion_lim = 1e6;
%[tempX,tempY] = ndgrid([CL.Bat,CL.Cap],[CL.Bat,CL.Cap]);
%all_diff_conn = [reshape(tempX,[size(tempX,1)*size(tempX,2),1]),reshape(tempY,[size(tempY,1)*size(tempY,2),1])];
all_diff_conn = nchoosek([CL.Bat],2);
linked_conn = {};
exhaust_num = nchoosek(size(all_diff_conn,1),converter_num);

% for i = 1:size(CL.Bat,2)  % use the averaged value of q for primary converter design
%     CL.Bat{i}.qlim = CL.Bat{i}.qlim_mu ;
% end
% CL = Flatten_battery_stat(CL); % use the flatened distributions of q for primary converter design
% CL = Update_cap_volt(CL);

if (exhaust_num<=trial_num)&&(exhaust_num<=comb_explosion_lim)
    all_diff_conn_choose_conv = nchoosek(1:size(all_diff_conn,1),converter_num);
    for j = 1:nchoosek(size(all_diff_conn,1),converter_num)
        linked_conn{end+1} =  all_diff_conn(all_diff_conn_choose_conv(j,:),:);  % if the search space is too large, do random sampling
    end    
    fprintf('%s%d','The search space',nchoosek(size(all_diff_conn,1),converter_num),'is complete');
    fprintf('\n');
else 
    for j = 1:trial_num
        linked_conn{end+1} = datasample(all_diff_conn,converter_num,'Replace',false);   % if the search space is too large, do random sampling
    end
    fprintf('%s%d%s%d','The search space ',nchoosek(size(all_diff_conn,1),converter_num),...
        ' is partially reached by ',trial_num);
    fprintf('\n');
end
%OptRes = repmat(struct('sol',[],'fval',[],'P_diff',[],'P_direct',[]),trial_num,1);
%% Solve the two-layer energy flow by parallel computing
for i = 1:CL.Stat.Bat_num
    Bat_Info(i) = CL.Bat{i}.qlim; 
end
parfor i = 1:min(trial_num,exhaust_num)
    P_diff_rating_mat_in = zeros(CL.Stat.Bat_num + CL.Stat.Bus_num, CL.Stat.Bat_num + CL.Stat.Bus_num);    % P_diff_mat_in connection should be reset in every cycle
    for j = 1:converter_num
         P_diff_rating_mat_in(linked_conn{i}{j,1}.ind, linked_conn{i}{j,2}.ind) = ENV.Avg_Conv.e_lim_vec(j);
    end
    OptRes(i) = func_dc_energyflow_two_layers(P_diff_rating_mat_in + P_diff_rating_mat_in', Bat_Info, CL, ENV);
end
%% Pick the Optimal Data
neg_max_q_out = [];
for i = 1:min(trial_num,exhaust_num)
    neg_max_q_out = [neg_max_q_out,OptRes(i).fval];
end
OutRes_energy_collection = {};
for i = 1:min(trial_num,exhaust_num)
    if (OptRes(i).fval == min(neg_max_q_out))
        
        OutRes_energy_collection{end+1}.P_diff_avg_mat_logic = OptRes(i).P_diff_avg;
        OutRes_energy_collection{end}.P_diff_avg_mat_rating = OptRes(i).Conv_energy_rating_partition_mat;
        OutRes_energy_collection{end}.Avg_energy_flow = OptRes(i).Avg_energy_flow;
        OutRes_energy_collection{end}.Var_energy_flow = OptRes(i).Var_energy_flow;
        OutRes_energy_collection{end}.Direct_energy_flow = OptRes(i).Direct_energy_flow;
 %       OutRes_energy_collection.Total_energy_flow = merge_multilayer_flow(OutRes_energy_collection.Avg_energy_flow,OutRes_energy_collection.Var_energy_flow,CL);
        
        OutRes_energy_collection{end}.Var_energy_process = OptRes(i).Var_energy_process;
        OutRes_energy_collection{end}.Total_energy_process = OptRes(i).Total_energy_process;
        OutRes_energy_collection{end}.Avg_energy_process = OptRes(i).Avg_energy_process;
        
        OutRes_energy_collection{end}.Maximum_output_energy = abs(OptRes(i).Maximum_output_energy);
 %       OutRes_energy_collection.Var_conv_energy_rating_vec = OptRes(i).Var_conv_energy_rating;
 %       OutRes_energy_collection.Avg_conv_energy_rating_vec = ENV.Avg_Conv.e_lim;
        
        OutRes_energy_collection{end}.Bat_energy_rating = OptRes(i).Bat_energy_rating;
        OutRes_energy_collection{end}.Bat_uratio = abs(OptRes(i).Maximum_output_energy)/sum(OptRes(i).Bat_energy_rating);
        OutRes_energy_collection{end}.CL_info = OptRes(i).Bat_energy_rating;
%        OptAvg.Conv_energy_rate_partition = func_rate_partition(OptRes(i).Conv_energy_rate,ENV.Avg_Conv.partition);         % sort the cnverter rating in ascending orders
%       fprintf('%s%d\n','Battery 1 to Battery 2: ',OptRes(i).sol((CL.Bat{2}.ind-1)*Node_num+CL.Bat{1}.ind));
%       fprintf('%s%d\n','Battery 2 to Battery 3: ',OptRes(i).sol((CL.Bat{3}.ind-1)*Node_num+CL.Bat{2}.ind));
%       fprintf('%s%d\n','Battery 3 to Capacitor 1: ',OptRes(i).sol((CL.Cap{1}.ind-1)*Node_num+CL.Bat{3}.ind));
%         fprintf('%s%d\n','Battery 1 to Bus 1: ',OptRes(i).sol((CL.Bus{1}.ind-1)*Node_num+CL.Bat{1}.ind));
%         fprintf('%s%d\n','Battery 2 to Bus 1: ',OptRes(i).sol((CL.Bus{1}.ind-1)*Node_num+CL.Bat{2}.ind));
%         fprintf('%s%d\n','Battery 3 to Bus 1: ',OptRes(i).sol((CL.Bus{1}.ind-1)*Node_num+CL.Bat{3}.ind));
%         fprintf('%s%d\n','Capacitor 1 to Bus 1: ',OptRes(i).sol((CL.Bus{1}.ind-1)*Node_num+CL.Cap{1}.ind));
    end
end

if(size(OutRes_energy_collection,2) > 1)
    fprintf('Multiple Optimal Points, use lexicographic optimization\n');
    % optimize for the processed energy
    for i = 1:size(OutRes_energy_collection,2)
        processed_energy(i) = sum(sum(OutRes_energy_collection{i}.Total_energy_process))/2;
    end
    OutRes_energy_processed_optimized = [];
    for i = 1:size(OutRes_energy_collection,2)
        if((sum(sum(OutRes_energy_collection{i}.Total_energy_process))/2) == min(processed_energy))
           OutRes_energy_processed_optimized{end+1} = OutRes_energy_collection{i};
        end
    end
     % if still multiple optimal, select by hands...
    if(size(OutRes_energy_processed_optimized,2) == 1)
        OutRes = OutRes_energy_processed_optimized{1};
    else
        save('Design_energy_collection', 'OutRes_energy_processed_optimized');
        fprintf('Still have multiple Optimal Points after lexicographic optimization\n');
        fprintf('Open Design_energy_collection.m \n');
        u_in = input('Assign xth entry of OutRes_energy_processed_optimized to OutRes, Type a number: \n');
        OutRes = OutRes_energy_processed_optimized{u_in};
    end
else
    OutRes = OutRes_energy_collection{1};
end

% function y = Flatten_battery_stat(CL)
%     temp = norminv(linspace(1e-4,1-1e-4,CL.Stat.Bat_num+1),CL.Bat{1}.qlim_mu, CL.Bat{1}.qlim_var);   %Normal inverse cumulative distribution function
%     for k = 1:CL.Stat.Bat_num
%        x_int = linspace(temp(k),temp(k+1),100);
%        y_int = normpdf(x_int,CL.Bat{1}.qlim_mu, CL.Bat{1}.qlim_var).*x_int;
%        CL.Bat{k}.qlim = CL.Stat.Bat_num*trapz(x_int,y_int);
%     end
%     y = CL;
% end
% 
% function y = Update_cap_volt(CL)
%    y = CL;
% end

end

