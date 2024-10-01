function OutRes = func_two_converter_power_design(CL,ENV)
%Input the diffferetial power delivery connection
trial_num = ENV.Avg_Conv.trial_num;
converter_num = ENV.Avg_Conv.Num;

comb_explosion_lim = 1e6;
%[tempX,tempY] = ndgrid([CL.Bat,CL.Cap],[CL.Bat,CL.Cap]);
%all_conn = [reshape(tempX,[size(tempX,1)*size(tempX,2),1]),reshape(tempY,[size(tempY,1)*size(tempY,2),1])];
all_diff_conn = nchoosek([CL.Bat],2);
linked_conn = {};
exhaust_num = nchoosek(size(all_diff_conn,1),converter_num);

% for i = 1:size(CL.Bat,2)  % use the averaged value of q for primary converter design
%     CL.Bat{i}.curlim = CL.Bat{i}.curlim_mu ;
% end
% CL = Flatten_battery_stat(CL); % use the flatened distributions of q for primary converter design
% CL = Update_cap_volt(CL);

if (exhaust_num<=trial_num)&&(exhaust_num<=comb_explosion_lim)
    all_conn_choose_conv = nchoosek(1:size(all_diff_conn,1),converter_num);
    for j = 1:nchoosek(size(all_diff_conn,1),converter_num)
        linked_conn{end+1} =  all_diff_conn(all_conn_choose_conv(j,:),:);  % if the search space is too large, do random sampling
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
%% Solve the two-layer power flow by parallel computing
for i = 1:CL.Stat.Bat_num
    Bat_Info(i) = CL.Bat{i}.curlim; 
end
parfor i = 1:min(trial_num,exhaust_num)
    P_diff_rating_mat_in = zeros(CL.Stat.Bat_num + CL.Stat.Bus_num, CL.Stat.Bat_num + CL.Stat.Bus_num);    % P_diff_mat_in connection should be reset in every cycle
    for j = 1:converter_num
         P_diff_rating_mat_in(linked_conn{i}{j,1}.ind, linked_conn{i}{j,2}.ind) = ENV.Avg_Conv.p_lim_vec(j);
    end
    OptRes(i) = func_dc_powerflow_two_layers(P_diff_rating_mat_in + P_diff_rating_mat_in', Bat_Info, CL, ENV);
end
%% Pick the Optimal Data
neg_max_q_out = [];
for i = 1:min(trial_num,exhaust_num)
    neg_max_q_out = [neg_max_q_out,OptRes(i).fval];
end

OutRes_power_collection = {};

for i = 1:min(trial_num,exhaust_num)
    if (OptRes(i).fval == min(neg_max_q_out))
        
        OutRes_power_collection{end+1}.P_diff_avg_mat_logic = OptRes(i).P_diff_avg;
        OutRes_power_collection{end}.P_diff_avg_mat_rating = OptRes(i).Conv_power_rating_partition_mat;
        OutRes_power_collection{end}.Avg_power_flow = OptRes(i).Avg_power_flow;
        OutRes_power_collection{end}.Var_power_flow = OptRes(i).Var_power_flow;
        OutRes_power_collection{end}.Direct_power_flow = OptRes(i).Direct_power_flow;
 %       OutRes_power_collection.Total_power_flow = merge_multilayer_flow(OutRes_power_collection.Avg_power_flow,OutRes_power_collection.Var_power_flow,CL);
        
        OutRes_power_collection{end}.Var_power_process = OptRes(i).Var_power_process;
        OutRes_power_collection{end}.Total_power_process = OptRes(i).Total_power_process;
        OutRes_power_collection{end}.Avg_power_process = OptRes(i).Avg_power_process;
        
        OutRes_power_collection{end}.Maximum_output_power = abs(OptRes(i).Maximum_output_power);
        
        OutRes_power_collection{end}.Bat_power_rating = OptRes(i).Bat_power_rating;
        OutRes_power_collection{end}.Bat_uratio = abs(OptRes(i).Maximum_output_power)/sum(OptRes(i).Bat_power_rating);
        OutRes_power_collection{end}.CL_info = OptRes(i).Bat_power_rating;
%        OutRes.Conv_power_rate_partition = func_rate_partition(OptRes(i).Conv_power_rate,ENV.Avg_Conv.partition);         % sort the cnverter rating in ascending orders
%       fprintf('%s%d\n','Battery 1 to Battery 2: ',OptRes(i).sol((CL.Bat{2}.ind-1)*Node_num+CL.Bat{1}.ind));
%       fprintf('%s%d\n','Battery 2 to Battery 3: ',OptRes(i).sol((CL.Bat{3}.ind-1)*Node_num+CL.Bat{2}.ind));
%       fprintf('%s%d\n','Battery 3 to Capacitor 1: ',OptRes(i).sol((CL.Cap{1}.ind-1)*Node_num+CL.Bat{3}.ind));
%         fprintf('%s%d\n','Battery 1 to Bus 1: ',OptRes(i).sol((CL.Bus{1}.ind-1)*Node_num+CL.Bat{1}.ind));
%         fprintf('%s%d\n','Battery 2 to Bus 1: ',OptRes(i).sol((CL.Bus{1}.ind-1)*Node_num+CL.Bat{2}.ind));
%         fprintf('%s%d\n','Battery 3 to Bus 1: ',OptRes(i).sol((CL.Bus{1}.ind-1)*Node_num+CL.Bat{3}.ind));
%         fprintf('%s%d\n','Capacitor 1 to Bus 1: ',OptRes(i).sol((CL.Bus{1}.ind-1)*Node_num+CL.Cap{1}.ind));
    end
end

if(size(OutRes_power_collection,2) > 1)
    fprintf('Multiple Optimal Points, use lexicographic optimization\n');
    % optimize for the processed power
    for i = 1:size(OutRes_power_collection,2)
        processed_power(i) = sum(sum(OutRes_power_collection{i}.Total_power_process))/2;
    end
    OutRes_power_processed_optimized = [];
    for i = 1:size(OutRes_power_collection,2)
        if((sum(sum(OutRes_power_collection{i}.Total_power_process))/2) == min(processed_power))
           OutRes_power_processed_optimized{end+1} = OutRes_power_collection{i};
        end
    end
     % if still multiple optimal, select by hands...
    if(size(OutRes_power_processed_optimized,2) == 1)
        OutRes = OutRes_power_processed_optimized{1};
    else
        save('Design_power_collection', 'OutRes_power_processed_optimized');
        fprintf('Still have multiple Optimal Points after lexicographic optimization\n');
        fprintf('Open Design_power_collection.m \n');
        u_in = input('Assign xth entry of OutRes_power_processed_optimized to OutRes, Type a number: \n');
        OutRes = OutRes_power_processed_optimized{u_in};
    end
else
    OutRes = OutRes_power_collection{1};
end


end

