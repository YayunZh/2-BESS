function OptRes = func_dc_powerflow_fpp(P_diff_rating_mat_in, Bat_info, CL, ENV)

% ENV.Avg_Conv.p_lim_mat = P_diff_mat_in;
Conv_power_rating_partition_mat = zeros(CL.Stat.Bat_num + CL.Stat.Bus_num, CL.Stat.Bat_num + CL.Stat.Bus_num);

CL = Update_Bat_stat(Bat_info, CL);
Bat_power_rating = zeros(CL.Stat.Bat_num,1);
for i = 1:CL.Stat.Bat_num
    Bat_power_rating(i) = CL.Bat{i}.volt*CL.Bat{i}.curlim;
end

eta = 1; % assuming the converter efficiency is 100% 
Node_num = CL.Stat.Bat_num + CL.Stat.Bus_num;

P_diff_avg_mat = (P_diff_rating_mat_in>0) | (P_diff_rating_mat_in>0)';    % the interconnection matrix must be symmetric
P_direct_mat = zeros(size(CL.Conn.direct));           % no direct power dilivery   
P_fpp_mat = CL.Conn.direct | CL.Conn.direct';
%Con_direct_num = sum(sum(P_direct_mat))/2;           % count the number of direct connections
%Con_diff_num = sum(sum(P_diff_avg_mat))/2;           % count the number of differential power delivery connections
P_diff_var_mat = CL.Conn.diff_var | CL.Conn.diff_var';  % the interconnection matrix for variance converter must be symmetric
Var_conv_num = ENV.Var_Conv.Num; % number of the var converters
%lambda = ENV.Constraint.lambda;  % ratio of processed power/output power
P_fpp_lim = ENV.Fpp_Conv.p_lim_singlevar;  % power limits of var converters
P_var_lim = ENV.Var_Conv.p_lim_singlevar;  % power limits of direc-power-processing converters
%% Construct the optimization problems and solve the linear programming problems 
P = optimvar('P',Node_num, Node_num);    % dimensional power transfer variables, primary layer
P_var = optimvar('P_var',Node_num, Node_num);    % 3*3 dimensional power transfer variables, secondary layer
% Xu = optimvar('Xu',Con_diff_num,'LowerBound',zeros(Con_diff_num,1));  % 3*1 dimensional non negative slack variables
bus_curlim = optimvar('bus_curlim',size(CL.Bus,2),1);     % bus current capability
bus_volt = [];
for i = 1:size(CL.Bus,2)
    CL.Bus{i}.curlim = bus_curlim(i);
    CL.Bus{i}.volt = 10;
%     for  j = 1:size(CL.Bat,2)
%        CL.Bus{i}.volt = CL.Bus{i}.volt + CL.Bat{j}.volt;
%     end
    bus_volt(i) = CL.Bus{i}.volt;
end
prob = optimproblem('Objective',sum(bus_curlim.*bus_volt),'ObjectiveSense','min');  % minimize the summations
% Construct the constraint for the slack power process ratio 
%power_ratio_con = optimconstr(1);
% power_ratio_con = Xu'*ones(Con_diff_num,1) <= -ENV.Constraint.lambda*sum(bus_curlim.*bus_volt); 
%prob.Constraints.power_ratio_con = power_ratio_con;

% Construct the constraint for direct power processing links  
% note that there is not any direct power processing links 
topo_direct_con = optimconstr(0);
% for i = 1:size(CL.Bat,2)
%     for j = 1:size(CL.Bus,2)
%         if (P_direct_mat(CL.Bat{i}.ind,CL.Bus{j}.ind) == 1)
%             topo_direct_con(end+1) = P(CL.Bat{i}.ind,CL.Bus{j}.ind) == - CL.Bat{i}.volt*CL.Bus{j}.curlim;
%         end
%         if (P_direct_mat(CL.Bus{j}.ind, CL.Bat{i}.ind) == 1)
%             topo_direct_con(end+1) = P(CL.Bus{j}.ind,CL.Bat{i}.ind) == CL.Bat{i}.volt*CL.Bus{j}.curlim;
%         end
%     end
% end
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

% Construct the constriant for full power processing links
topo_fpp_con = optimconstr(0);
for i = 1:Node_num
    for j = i:Node_num
        if (P_fpp_mat(i,j) == 1)
            topo_fpp_con(end+1) = P(i,j) == -P(j,i);
        else
            topo_fpp_con(end+1) = P(i,j) == 0;
            topo_fpp_con(end+1) = P(j,i) == 0; 
        end
    end
end
 prob.Constraints.topo_fpp_con = topo_fpp_con;

% Construct the constriant for differential power processing links, primary
% (average converters)
topo_diff_avg_con = optimconstr(0);
for i = 1:Node_num
    for j = i:Node_num
        if (P_fpp_mat(i,j) == 0)
            if (P_diff_avg_mat(i,j) == 1)
                topo_diff_avg_con(end+1) = P(i,j) == -P(j,i);
            else
                topo_diff_avg_con(end+1) = P(i,j) == 0;
                topo_diff_avg_con(end+1) = P(j,i) == 0; 
            end
        end
    end
end
 prob.Constraints.topo_diff_avg_con = topo_diff_avg_con;
 
% Construct the constriant for differential power processing links,
% secondary (var converters)
topo_diff_var_con = optimconstr(0);
for i = 1:Node_num
    for j = i:Node_num
        if (P_diff_var_mat(i,j) == 1)
            topo_diff_var_con(end+1) = P_var(i,j) == -P_var(j,i);
        else
            topo_diff_var_con(end+1) = P_var(i,j) == 0;
            topo_diff_var_con(end+1) = P_var(j,i) == 0; 
        end
    end
end
 prob.Constraints.topo_diff_var_con = topo_diff_var_con;

% Construct the primary/average converter power limits
avg_conv_con = optimconstr(0);
temp_ct = 1;
for i = 1:Node_num
    for j = i:Node_num
        if (P_diff_avg_mat(i,j) == 1)
            avg_conv_con(end+1) = P(i,j) <=  P_diff_rating_mat_in(i,j);
            avg_conv_con(end+1) = -P_diff_rating_mat_in(i,j) <= P(i,j);
%            Conv_power_rating_partition_mat(i,j) = ENV.Avg_Conv.p_lim(temp_ct);
        end
    end
end
prob.Constraints.avg_conv_con = avg_conv_con; 

% Construct the full-power-processing converter power limits
fpp_conv_con = optimconstr(0);
for i = 1:Node_num
    for j = i:Node_num
        if (P_fpp_mat(i,j) == 1)
            fpp_conv_con(end+1) = P(i,j) <= P_fpp_lim;
            fpp_conv_con(end+1) = -P_fpp_lim <= P(i,j);
%            Conv_power_rating_partition_mat(i,j) = ENV.Avg_Conv.p_lim(temp_ct);
        end
    end
end
prob.Constraints.fpp_conv_con = fpp_conv_con;

% Construct the primary/average converter power limits
avg_conv_con = optimconstr(0);
temp_ct = 1;
for i = 1:Node_num
    for j = i:Node_num
        if (P_diff_avg_mat(i,j) == 1)
            avg_conv_con(end+1) = P(i,j) <=  P_diff_rating_mat_in(i,j);
            avg_conv_con(end+1) = -P_diff_rating_mat_in(i,j) <= P(i,j);
%            Conv_power_rating_partition_mat(i,j) = ENV.Avg_Conv.p_lim(temp_ct);
        end
    end
end
prob.Constraints.avg_conv_con = avg_conv_con;

% Construct the power limit for secondary/variance power converters
var_conv_con = optimconstr(0);
for i = 1:Node_num
    for j = i:Node_num
        if (P_diff_var_mat(i,j) == 1)
            var_conv_con(end+1) = P_var(i,j) <= P_var_lim;
            var_conv_con(end+1) = -P_var_lim <=  P_var(i,j);
        end
    end
end
prob.Constraints.var_conv_con = var_conv_con;

% Construct the constriant for power conservation
power_bus_con = optimconstr(0);
for i = 1:size(CL.Bus,2)
    power_bus_con(end+1) = sum(P(CL.Bus{i}.ind,:)) == CL.Bus{i}.volt*CL.Bus{i}.curlim;
end
prob.Constraints.power_bus_con = power_bus_con;

% power_cap_con = optimconstr(0);
% for i = 1:size(CL.Cap,2)
%     power_cap_con(end+1) = (sum(P(CL.Cap{i}.ind,:)) + sum(P_var(CL.Cap{i}.ind,:))) == 0;
% end
% prob.Constraints.power_cap_con = power_cap_con;

power_bat_con = optimconstr(0);
for i = 1:size(CL.Bat,2)
    power_bat_con(end+1) = CL.Bat{i}.volt*CL.Bat{i}.curlim >= sum(P(CL.Bat{i}.ind,:));
    power_bat_con(end+1) = sum(P(CL.Bat{i}.ind,:)) >= -CL.Bat{i}.volt*CL.Bat{i}.curlim;
end
prob.Constraints.power_bat_con = power_bat_con;   


problem = prob2struct(prob);
[sol,fval,exitflag,output] = linprog(problem);

if(exitflag == 1)
    OptRes.sol = sol;
    OptRes.fval = fval;
    OptRes.P_diff_avg = P_diff_avg_mat;
    OptRes.P_diff_var = P_diff_var_mat;
    OptRes.P_direct = P_direct_mat;
    OptRes.P_full_power_processing = P_fpp_mat;
    OptRes.Conv_power_rating_partition_mat = P_diff_rating_mat_in;
    
    Total_power_flow = reshape(OptRes.sol(1:Node_num*Node_num),[Node_num,Node_num]);
    [OptRes.Avg_power_flow, OptRes.Fpp_power_flow] = func_power_flow_decomp(Total_power_flow, CL);
    OptRes.Var_power_flow = reshape(OptRes.sol(Node_num*Node_num+1:2*Node_num*Node_num),[Node_num,Node_num]);
    
    OptRes.Var_power_process = sum(sum(abs(OptRes.Var_power_flow)))/2;
    OptRes.Avg_power_process = sum(sum(abs(OptRes.Avg_power_flow)))/2;
    OptRes.Fpp_power_process = sum(sum(abs(OptRes.Fpp_power_flow)))/2;
    OptRes.Total_power_process = OptRes.Fpp_power_process + OptRes.Var_power_process + OptRes.Avg_power_process;
    
    OptRes.ComponentInfo = CL;

    OptRes.Maximum_output_power = fval;
    for i = 1: CL.Stat.Bat_num
        OptRes.BatOutputpower(i) = sum((OptRes.Var_power_flow(CL.Bat{i}.ind,:))) + ...
         sum((OptRes.Avg_power_flow(CL.Bat{i}.ind,:))) + sum((OptRes.Fpp_power_flow(CL.Bat{i}.ind,:)));
    end
    OptRes.Var_conv_power_rating = ones(1,Var_conv_num)*P_var_lim;
    OptRes.Bat_power_rating = Bat_power_rating;
    
    if (sol(end) == sol(2*Node_num*Node_num + 1))
        OptRes.Output_q_max = sol(2*Node_num*Node_num + 1:end);
    else
        error('Error in interpreting results in dc power flow optimization of two lay converters\n');
    end
else
    error('Error in dc power flow optimization on two-layer converters\n');
end

function y = Update_Bat_stat(Bat_Info,CL)
    for l = 1:CL.Stat.Bat_num
        %CL.Bat{l}.curlim = mvnrnd([CL.Bat{l}.curlim_mu,CL.Bat{l}.volt_mu], diag([CL.Bat{l}.curlim_var,CL.Bat{l}.volt_var]),1);
         CL.Bat{l}.curlim = Bat_Info(l);
    end
    CL.Cap{1}.volt = CL.Bus{1}.volt - CL.Bat{1}.volt - CL.Bat{2}.volt - CL.Bat{3}.volt...
    - CL.Bat{4}.volt - CL.Bat{5}.volt - CL.Bat{6}.volt...
    - CL.Bat{7}.volt - CL.Bat{8}.volt - CL.Bat{9}.volt;     % Decide the capactior votlages
    y = CL;
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

