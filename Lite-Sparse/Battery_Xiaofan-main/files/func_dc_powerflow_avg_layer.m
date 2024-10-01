function OptRes = func_ac_powerflow_avg_layer(P_diff_mat_in, CL, ENV)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
P_diff_mat = (P_diff_mat_in>0) | (P_diff_mat_in>0)';           % the interconnection matrix must be symmetric
Con_diff_num = sum(sum(P_diff_mat))/2;      % count the number of differential power delivery connections
eta = 1; % assuming the converter efficiency is 100% 
Conv_power_rating_partition_mat = zeros(CL.Stat.Bat_num,CL.Stat.Bat_num);
Node_num = CL.Stat.Bat_num;
%% Construct the optimization problems and solve the linear programming problems 
I_B = zeros(n_bat,n_t); % battery charge
I_C_s = zeros(n_bat-1,n_t); % sparse converter charge
I_C_d = zeros(n_bat-1,n_t); % dense converter charge
I_L = zeros(n_bat, n_t); % load charge

P = optimvar('P',Node_num, Node_num);    % dimensional power transfer variables
%Xu = optimvar('Xu',Con_diff_num,'LowerBound',zeros(Con_diff_num,1));  % 3*1 dimensional non negative slack variables
bus_curlim = optimvar('bus_curlim',size(CL.Bus,2),1);     % bus current capability
bus_volt = [];
for i = 1:size(CL.Bus,2)
    CL.Bus{i}.curlim = bus_curlim(i);
    CL.Bus{i}.volt = 0;
    for  j = 1:size(CL.Bat,2)
       CL.Bus{i}.volt = CL.Bus{i}.volt + CL.Bat{j}.volt;
    end
    bus_volt(i) = CL.Bus{i}.volt;
end
prob = optimproblem('Objective',sum(bus_curlim.*bus_volt),'ObjectiveSense','min');  % minimize the summations
% Construct the constraint for the slack power process ratio 
% power_ratio_con = optimconstr(1);
% power_ratio_con = Xu'*ones(Con_diff_num,1) <= -ENV.Constraint.lambda*sum(bus_curlim.*bus_volt); 
% prob.Constraints.power_ratio_con = power_ratio_con;
% Construct the constraint for direct power processing links
topo_direct_con = optimconstr(0);
for i = 1:size(CL.Bat,2)
    for j = 1:size(CL.Bus,2)
        if (P_direct_mat(CL.Bat{i}.ind,CL.Bus{j}.ind) == 1)
            topo_direct_con(end+1) = P(CL.Bat{i}.ind,CL.Bus{j}.ind) == - CL.Bat{i}.volt*CL.Bus{j}.curlim;
        end
        if (P_direct_mat(CL.Bus{j}.ind, CL.Bat{i}.ind) == 1)
            topo_direct_con(end+1) = P(CL.Bus{j}.ind,CL.Bat{i}.ind) == CL.Bat{i}.volt*CL.Bus{j}.curlim;
        end
    end
end
% for i = 1:size(CL.Cap,2)
%     for j = 1:size(CL.Bus,2)
%         if (P_direct_mat(CL.Cap{i}.ind,CL.Bus{j}.ind) == 1)
%             topo_direct_con(end+1) = P(CL.Cap{i}.ind, CL.Bus{j}.ind) == - CL.Cap{i}.volt*CL.Bus{j}.curlim;
%         end
%         if (P_direct_mat(CL.Bus{j}.ind, CL.Cap{i}.ind) == 1)
%             topo_direct_con(end+1) = P(CL.Bus{j}.ind, CL.Cap{i}.ind) == CL.Cap{i}.volt*CL.Bus{j}.curlim;
%         end 
%     end
% end
prob.Constraints.topo_direct_con = topo_direct_con;

% Construct the constriant for differential power processing links
topo_diff_con = optimconstr(0);
for i = 1:Node_num
    for j = i:Node_num
        if (P_direct_mat(i,j) == 0)
            if (P_diff_mat(i,j) == 1)
                topo_diff_con(end+1) = P(i,j) == -P(j,i);
            else
                topo_diff_con(end+1) = P(i,j) == 0;
                topo_diff_con(end+1) = P(j,i) == 0; 
            end
        end
    end
end
 prob.Constraints.topo_diff_con = topo_diff_con;
 % Construct the slack var
%  slack_con = optimconstr(0);
%  ct = 1;
% for i = 1:Node_num
%     for j = i:Node_num
%         if (P_diff_mat(i,j) == 1)
%             slack_con(end+1) = P(i,j) <= Xu(ct);
%             slack_con(end+1) = -Xu(ct) <= P(i,j);
%             ct = ct + 1;
%         end
%     end
% end
% if (ct ~= Con_diff_num + 1)
%     fprintf('slack constraint error!\n')
% else
%     prob.Constraints.slack_con = slack_con;
% end 
% Construct the Average Converter Power Limit
avg_conv_con = optimconstr(0);
temp_ct = 1;
for i = 1:Node_num
    for j = i:Node_num
        if (P_diff_mat(i,j) == 1)
%             avg_conv_con(end+1) = P(i,j) <= ENV.Avg_Conv.e_lim(temp_ct);
%             avg_conv_con(end+1) = -ENV.Avg_Conv.e_lim(temp_ct) <= P(i,j);
%             Conv_power_rating_partition_mat(i,j) = ENV.Avg_Conv.e_lim(temp_ct);
            avg_conv_con(end+1) = P(i,j) <= P_diff_mat_in(i,j);
            avg_conv_con(end+1) = -P_diff_mat_in(i,j) <= P(i,j);
            Conv_power_rating_partition_mat(i,j) = P_diff_mat_in(i,j);
            temp_ct = temp_ct + 1;
        end
    end
end
prob.Constraints.avg_conv_con = avg_conv_con;

% Construct the constriant for power conservation
power_bus_con = optimconstr(0);
for i = 1:size(CL.Bus,2)
    power_bus_con(end+1) = sum(P(CL.Bus{i}.ind,:))== CL.Bus{i}.volt*CL.Bus{i}.curlim;
end
prob.Constraints.power_bus_con = power_bus_con;

% power_cap_con = optimconstr(0);
% for i = 1:size(CL.Cap,2)
%     power_cap_con(end+1) = sum(P(CL.Cap{i}.ind,:))== 0;
% end
% prob.Constraints.power_cap_con = power_cap_con;

power_bat_con = optimconstr(0);
for i = 1:size(CL.Bat,2)
    power_bat_con(end+1) = CL.Bat{i}.volt*CL.Bat{i}.curlim >= sum(P(CL.Bat{i}.ind,:));
    power_bat_con(end+1) = sum(P(CL.Bat{i}.ind,:))>= -CL.Bat{i}.volt*CL.Bat{i}.curlim;
end
prob.Constraints.power_bat_con = power_bat_con;

problem = prob2struct(prob);
[sol,fval,exitflag,output] = linprog(problem);

if(exitflag == 1)
    OptRes.sol = sol;
    OptRes.bus_plim = fval;
    OptRes.Conv_power_rating_partition_mat = Conv_power_rating_partition_mat + Conv_power_rating_partition_mat';
%    OptRes.P_diff = P_diff_mat;
%    OptRes.P_direct = P_direct_mat;
    Total_power_flow = reshape(OptRes.sol(1:Node_num*Node_num),[Node_num,Node_num]);
    [OptRes.Avg_power_flow, OptRes.Direct_power_flow] = func_power_flow_decomp(Total_power_flow,CL);
%    OptRes.diff_mat_avg = P_diff_mat;
%    OptRes.Conv_power_processed = abs(sol((Node_num*Node_num+1):(Node_num*Node_num+Con_diff_num)));
    if (sol(end) == sol(Node_num*Node_num+1))
        OptRes.Output_cur_max = sol(Node_num*Node_num+1:end);
    else
        error('Error in interpreting results in dc power flow optimization of average converters\n');
    end
else
    error('Error in optimizing dc power flow of the averaging converters\n');
end

function [Indirect_power_flow,Direct_power_flow] = func_power_flow_decomp(Total_power_flow,CL)
    Indirect_power_flow = zeros(size(Total_power_flow,1),size(Total_power_flow,2));
    Direct_power_flow = zeros(size(Total_power_flow,1),size(Total_power_flow,2));
    for ic = 1:size(Total_power_flow,1)
        for jc = 1:size(Total_power_flow,2)
           if ((ic < CL.Bus{1}.ind)&&(jc < CL.Bus{1}.ind))
               Indirect_power_flow(ic,jc) = Total_power_flow(ic,jc);
           else
               Direct_power_flow(ic,jc) = Total_power_flow(ic,jc);
           end
        end
    end 
end

end

