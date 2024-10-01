function OptAvg = func_average_converter_energy_design(CL,ENV)
%Input the differential energy delivery connection
trial_num = ENV.Avg_Conv.trial_num;
converter_num = ENV.Avg_Conv.Num;
P_direct_mat_in = CL.Conn.direct;
Node_num = CL.Stat.Bat_num + CL.Stat.Bus_num;

comb_explosion_lim = 1e6;
%[tempX,tempY] = ndgrid([CL.Bat,CL.Cap],[CL.Bat,CL.Cap]);
%all_diff_conn = [reshape(tempX,[size(tempX,1)*size(tempX,2),1]),reshape(tempY,[size(tempY,1)*size(tempY,2),1])];
all_diff_conn = nchoosek([CL.Bat],2);% Take any combination of two of them, 36*2 cell
linked_conn = {};
exhaust_num = nchoosek(size(all_diff_conn,1),converter_num); % the possible place of avg_conv

for i = 1:size(CL.Bat,2)  % use the averaged value of q for primary converter design
    CL.Bat{i}.qlim = CL.Bat{i}.qlim_mu ;
end

% battery creation
CL = Flatten_battery_stat(CL); % use the flatened distributions of q for primary converter design
CL = Updat_cap_volt(CL);

Bat_energy_rating = zeros(CL.Stat.Bat_num,1);
for i = 1:CL.Stat.Bat_num 
    Bat_energy_rating(i) = CL.Bat{i}.volt*CL.Bat{i}.qlim;
end

% Iterate over all possibilities
% Or random sampling if the number of combinations is too large
if (exhaust_num<=trial_num)&&(exhaust_num<=comb_explosion_lim)
    all_diff_conn_choose_conv = nchoosek(1:size(all_diff_conn,1),converter_num);% 7140, positions of 3 conv
    for j = 1:nchoosek(size(all_diff_conn,1),converter_num)% 7140
        % every element: 3*2 cell, only between batteries
        linked_conn{end+1} =  all_diff_conn(all_diff_conn_choose_conv(j,:),:);  
    end    
    fprintf('%s%d','The search space ',nchoosek(size(all_diff_conn,1),converter_num),' is complete');
    fprintf('\n');
else 
    for j = 1:trial_num
        linked_conn{end+1} = datasample(all_diff_conn,converter_num,'Replace',false);   % if the search space is too large, do random sampling
    end
    fprintf('%s%d%s%d','The search space ',nchoosek(size(all_diff_conn,1),converter_num),...
        ' is partially reached by ',trial_num);
    fprintf('\n');
end
%% Find the optimal topology by parallel computing
%parfor i = 1:min(trial_num,exhaust_num) % parallel loop
for i = 1:min(trial_num,exhaust_num)
    P_diff_mat_in = zeros(CL.Stat.Bat_num + CL.Stat.Bus_num, CL.Stat.Bat_num + CL.Stat.Bus_num);    % P_diff_mat_in connection should be reset in every cycle
    for j = 1:converter_num
         P_diff_mat_in(linked_conn{i}{j,1}.ind, linked_conn{i}{j,2}.ind) = ENV.Avg_Conv.e_lim_vec(j); % Index of start and end points
    end
    % 3 conv setting done
    OptRes(i) = func_dc_energyflow_avg_layer(P_diff_mat_in + P_diff_mat_in', CL, ENV);
end
%% Pick the Optimal Data
neg_max_e_out = [];
for i = 1:min(trial_num,exhaust_num) % extract bus_elim
    neg_max_e_out = [neg_max_e_out,OptRes(i).bus_elim]; % sum(bus_qlim.*bus_volt)
end
OptAvg_energy_collection = {};
for i = 1:min(trial_num,exhaust_num)
    if (OptRes(i).bus_elim == min(neg_max_e_out)) % objective
        OptAvg_energy_collection{end+1}.output_energy = abs(OptRes(i).bus_elim);
        OptAvg_energy_collection{end}.Avg_energy_flow = OptRes(i).Avg_energy_flow; % bat - bat converter
        OptAvg_energy_collection{end}.Direct_energy_flow = OptRes(i).Direct_energy_flow; % bat - bus
        
        OptAvg_energy_collection{end}.diff_avg_mat = (OptRes(i).Avg_energy_flow ~= 0); % avg_converter connection
        OptAvg_energy_collection{end}.Conv_energy_rating_partition_mat = OptRes(i).Conv_energy_rating_partition_mat; %P_diff_mat_in
        OptAvg_energy_collection{end}.Avg_conv_energy_processed = abs(OptRes(i).Avg_energy_flow);

        OptAvg_energy_collection{end}.Bat_energy_rating = Bat_energy_rating;
        OptAvg_energy_collection{end}.Bat_uratio = abs(OptRes(i).bus_elim)/sum(Bat_energy_rating);
        OptAvg_energy_collection{end}.Avg_conv_uratio = sum(sum(abs(OptRes(i).Avg_energy_flow)))/sum(sum(OptRes(i).Conv_energy_rating_partition_mat));
%        OptAvg_energy_collection{end+1}.Conv_energy_rate_partition = func_rate_partition(OptRes(i).Conv_energy_rate,ENV.Avg_Conv.partition);         % sort the cnverter rating in ascending orders
      % fprintf('%s%d\n','Battery 1 to Battery 2: ',OptRes(i).sol((CL.Bat{2}.ind-1)*Node_num+CL.Bat{1}.ind));
      % fprintf('%s%d\n','Battery 2 to Battery 3: ',OptRes(i).sol((CL.Bat{3}.ind-1)*Node_num+CL.Bat{2}.ind));
      % fprintf('%s%d\n','Battery 3 to Capacitor 1: ',OptRes(i).sol((CL.Cap{1}.ind-1)*Node_num+CL.Bat{3}.ind));
      %   fprintf('%s%d\n','Battery 1 to Bus 1: ',OptRes(i).sol((CL.Bus{1}.ind-1)*Node_num+CL.Bat{1}.ind));
      %   fprintf('%s%d\n','Battery 2 to Bus 1: ',OptRes(i).sol((CL.Bus{1}.ind-1)*Node_num+CL.Bat{2}.ind));
      %   fprintf('%s%d\n','Battery 3 to Bus 1: ',OptRes(i).sol((CL.Bus{1}.ind-1)*Node_num+CL.Bat{3}.ind));
      %   fprintf('%s%d\n','Capacitor 1 to Bus 1: ',OptRes(i).sol((CL.Bus{1}.ind-1)*Node_num+CL.Cap{1}.ind));
    end
end

if(size(OptAvg_energy_collection,2) > 1)
    fprintf('Multiple Optimal Points, use lexicographic optimization\n');
    % optimize for the converter rating
    for i = 1:size(OptAvg_energy_collection,2)
        converter_rating(i) = max(max(OptAvg_energy_collection{i}.Avg_conv_energy_processed));
    end
    OptAvg_energy_rating_optimized = [];
    for i = 1:size(OptAvg_energy_collection,2)
        if(max(max(OptAvg_energy_collection{i}.Avg_conv_energy_processed)) == min(converter_rating))
           OptAvg_energy_rating_optimized{end+1} = OptAvg_energy_collection{i};
        end
    end
    % optimize for the processed energy
    for i = 1:size(OptAvg_energy_rating_optimized,2)
        processed_energy(i) = sum(sum(OptAvg_energy_rating_optimized{i}.Avg_conv_energy_processed))/2;
    end
    OptAvg_energy_processed_optimized = [];
    for i = 1:size(OptAvg_energy_rating_optimized,2)
        if((sum(sum(OptAvg_energy_rating_optimized{i}.Avg_conv_energy_processed))/2) == min(processed_energy))
           OptAvg_energy_processed_optimized{end+1} = OptAvg_energy_rating_optimized{i};
        end
    end
     % if still multiple optimal, select by hands...
    if(size(OptAvg_energy_processed_optimized,2)==1)
        OptAvg = OptAvg_energy_processed_optimized{1};
    else
        save('PreDesign_energy_collection', 'OptAvg_energy_processed_optimized');
        fprintf('Still have multiple Optimal Points after lexicographic optimization\n');
        fprintf('Open PreDesign_energy_collection.m \n');
        fprintf('The number of optimal points is: \n');
        size(OptAvg_energy_processed_optimized,2)
        u_in = input('Assign xth entry of OptAvg_energy_processed_optimized to OptAvg, Type a number: \n');
        OptAvg = OptAvg_energy_processed_optimized{u_in};
    end
else
    OptAvg = OptAvg_energy_collection{1};
end

function y = Flatten_battery_stat(CL)
    temp = norminv(linspace(1e-4,1-1e-4,CL.Stat.Bat_num+1),CL.Bat{1}.qlim_mu, CL.Bat{1}.qlim_var);   %Normal inverse cumulative distribution function
    for k = 1:CL.Stat.Bat_num
       x_int = linspace(temp(k),temp(k+1),100);
       y_int = normpdf(x_int,CL.Bat{1}.qlim_mu, CL.Bat{1}.qlim_var).*x_int;
       CL.Bat{k}.qlim = CL.Stat.Bat_num*trapz(x_int,y_int);
    end
    y = CL;
%     temp = norminv(linspace(1e-4,1-1e-4,CL.Stat.Bat_num+1),CL.Bat{1}.qlim_mu, CL.Bat{1}.qlim_var);   %Normal inverse cumulative distribution function
%     N_int = 100;
%     for k = 1:CL.Stat.Bat_num
%        x_int = linspace(temp(k),temp(k+1),N_int);
%        y_int = CL.Bat{1}.qcdf(x_int);
%        % we do integral by part, int xp(x)dx = delta xc(x) - int c(x)dx
%        CL.Bat{k}.qlim = ((x_int(N_int)*CL.Bat{1}.qcdf(x_int(N_int)) - ...
%            x_int(1)*CL.Bat{1}.qcdf(x_int(1)) - trapz(x_int, y_int)))/(1/CL.Stat.Bat_num);
%        bat_vec(k) = CL.Bat{k}.curlim;
%     end
%     y = CL;
end

function y = Update_cap_volt(CL)
   y = CL;
end

end

