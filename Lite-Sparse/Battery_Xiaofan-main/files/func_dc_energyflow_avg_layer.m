function OptRes = func_dc_energyflow_avg_layer(P_diff_mat_in, CL, ENV)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
% P_diff_mat_in: e_limit of 3 converters
Node_num = CL.Stat.Bat_num + CL.Stat.Bus_num;
% layer 1
P_diff_mat = (P_diff_mat_in>0) | (P_diff_mat_in>0)';  % the interconnection matrix must be symmetric
% B1-9 + C with Bus, 10 * 10
% direct should not be symmetric: P_direct_mat = CL.Conn.direct | CL.Conn.direct'
% prob. with end line, assume it as row bus
P_direct_mat = CL.Conn.direct | CL.Conn.direct';% the interconnection matrix must be symmetric
% P_direct_mat(end,:) = zeros(1,10);

% Con_direct_num calculation
Con_direct_num = sum(sum(P_direct_mat))/2;  % count the number of direct connections
Con_diff_num = sum(sum(P_diff_mat))/2;% 3 avg conv   % count the number of differential energy delivery connections
eta = 1; % assuming the converter efficiency is 100% 
Conv_energy_rating_partition_mat = zeros(CL.Stat.Bat_num + CL.Stat.Cap_num,CL.Stat.Bat_num + CL.Stat.Cap_num);
%% Construct the optimization problems and solve the linear programming problems 
P = optimvar('P',Node_num, Node_num);    % dimensional energy transfer variables
% Xu = optimvar('Xu',Con_diff_num,'LowerBound',zeros(Con_diff_num,1));  % 3*1 dimensional non negative slack variables
bus_qlim = optimvar('bus_qlim',size(CL.Bus,2),1);     % bus current capability
bus_volt = [];
% calculate CL.Bus{i}.volt, then give it to bus_volt
% bus_qlim :optivar; bus_volt: 9
for i = 1:size(CL.Bus,2)
    CL.Bus{i}.qlim = bus_qlim(i); %variable
    CL.Bus{i}.volt = 0;
    for  j = 1:size(CL.Bat,2)
        % v_bus = sum(v_b1-9)
       CL.Bus{i}.volt = CL.Bus{i}.volt + CL.Bat{j}.volt;
    end
    bus_volt(i) = CL.Bus{i}.volt;
end
% why minimum: because power output < 0
prob = optimproblem('Objective',sum(bus_qlim.*bus_volt),'ObjectiveSense','min');  % minimize the summations
% ? slack energy
% Construct the constraint for the slack energy process ratio 
% power_ratio_con = optimconstr(1);
% power_ratio_con = Xu'*ones(Con_diff_num,1) <= -ENV.Constraint.lambda*sum(bus_qlim.*bus_volt); 
% prob.Constraints.power_ratio_con = power_ratio_con;
% Construct the constraint for direct power processing links
topo_direct_con = optimconstr(0);
for i = 1:size(CL.Bat,2)% 9
    for j = 1:size(CL.Bus,2)% 1
        % sum = 0
        if (P_direct_mat(CL.Bat{i}.ind,CL.Bus{j}.ind) == 1)  % 1-9,10
            % constraints
            topo_direct_con(end+1) = P(CL.Bat{i}.ind,CL.Bus{j}.ind) == - CL.Bat{i}.volt*CL.Bus{j}.qlim;
        end
        if (P_direct_mat(CL.Bus{j}.ind, CL.Bat{i}.ind) == 1) % 10,1-9
            topo_direct_con(end+1) = P(CL.Bus{j}.ind,CL.Bat{i}.ind) == CL.Bat{i}.volt*CL.Bus{j}.qlim;
        end
    end
end
% for i = 1:size(CL.Cap,2)
%     for j = 1:size(CL.Bus,2)
%         if (P_direct_mat(CL.Cap{i}.ind,CL.Bus{j}.ind) == 1)
%             topo_direct_con(end+1) = P(CL.Cap{i}.ind, CL.Bus{j}.ind) == - CL.Cap{i}.volt*CL.Bus{j}.qlim;
%         end
%         if (P_direct_mat(CL.Bus{j}.ind, CL.Cap{i}.ind) == 1)
%             topo_direct_con(end+1) = P(CL.Bus{j}.ind, CL.Cap{i}.ind) == CL.Cap{i}.volt*CL.Bus{j}.qlim;
%         end 
%     end
% end
prob.Constraints.topo_direct_con = topo_direct_con;

% Construct the constriant for differential power processing links
topo_diff_con = optimconstr(0);
for i = 1:Node_num % 10
    for j = i:Node_num % 10
        if (P_direct_mat(i,j) == 0)
            if (P_diff_mat(i,j) == 1) % 3 avg conv
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
        if (P_diff_mat(i,j) == 1) % 3 avg conv
%             avg_conv_con(end+1) = P(i,j) <= ENV.Avg_Conv.e_lim(temp_ct);
%             avg_conv_con(end+1) = -ENV.Avg_Conv.e_lim(temp_ct) <= P(i,j);
%             Conv_energy_rating_partition_mat(i,j) = ENV.Avg_Conv.e_lim(temp_ct);
            % power constraints for each battery
            avg_conv_con(end+1) = P(i,j) <= P_diff_mat_in(i,j);
            avg_conv_con(end+1) = -P_diff_mat_in(i,j) <= P(i,j); % P(i,j) = P_bat
            Conv_energy_rating_partition_mat(i,j) = P_diff_mat_in(i,j);
            temp_ct = temp_ct + 1;
        end
    end
end
prob.Constraints.avg_conv_con = avg_conv_con;

% Construct the constriant for energy conservation
energy_bus_con = optimconstr(0);
for i = 1:size(CL.Bus,2)
    energy_bus_con(end+1) = sum(P(CL.Bus{i}.ind,:)) == CL.Bus{i}.volt*CL.Bus{i}.qlim;
end
prob.Constraints.energy_bus_con = energy_bus_con;

% energy_cap_con = optimconstr(0);
% for i = 1:size(CL.Cap,2)
%     energy_cap_con(end+1) = sum(P(CL.Cap{i}.ind,:))== 0;
% end
% prob.Constraints.energy_cap_con = energy_cap_con;

energy_bat_con = optimconstr(0);
for i = 1:size(CL.Bat,2)
    % only consider battery power output
    energy_bat_con(end+1) = CL.Bat{i}.volt*CL.Bat{i}.qlim >= sum(P(CL.Bat{i}.ind,:));
    energy_bat_con(end+1) = sum(P(CL.Bat{i}.ind,:))>= -CL.Bat{i}.volt*CL.Bat{i}.qlim;
end
prob.Constraints.energy_bat_con = energy_bat_con;

problem = prob2struct(prob);
[sol,fval,exitflag,output] = linprog(problem);

if(exitflag == 1)
    OptRes.sol = sol;
    OptRes.bus_elim = fval;
    OptRes.Conv_energy_rating_partition_mat = Conv_energy_rating_partition_mat + Conv_energy_rating_partition_mat';
%    OptRes.P_diff = P_diff_mat;
%    OptRes.P_direct = P_direct_mat;
    Total_energy_flow = reshape(OptRes.sol(1:Node_num*Node_num),[Node_num,Node_num]);
    [OptRes.Avg_energy_flow, OptRes.Direct_energy_flow] = func_energy_flow_decomp(Total_energy_flow,CL);
%    OptRes.diff_mat_avg = P_diff_mat;
%    OptRes.Conv_energy_processed = abs(sol((Node_num*Node_num+1):(Node_num*Node_num+Con_diff_num)));
    if (sol(end) == sol(Node_num*Node_num+1))
        OptRes.Output_q_max = sol(Node_num*Node_num+1:end);
    else
        error('Error in interpreting results in LP\n');
    end
else
    error('Error in designing averaging converters\n');
end

function [Indirect_energy_flow,Direct_energy_flow] = func_energy_flow_decomp(Total_energy_flow,CL)
    Indirect_energy_flow = zeros(size(Total_energy_flow,1),size(Total_energy_flow,2));
    Direct_energy_flow = zeros(size(Total_energy_flow,1),size(Total_energy_flow,2));
    for ic = 1:size(Total_energy_flow,1)% 10
        for jc = 1:size(Total_energy_flow,2)% 10
           if ((ic < CL.Bus{1}.ind)&&(jc < CL.Bus{1}.ind))% 1-9,1-9
               Indirect_energy_flow(ic,jc) = Total_energy_flow(ic,jc);
           else
               Direct_energy_flow(ic,jc) = Total_energy_flow(ic,jc);
           end
        end
    end 
end

end

