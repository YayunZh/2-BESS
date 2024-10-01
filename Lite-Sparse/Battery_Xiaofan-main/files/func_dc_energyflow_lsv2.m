function OptRes = func_dc_energyflow_lsv2(k2,Bat_info, CL, ENV)
% P_diff_rating_mat_in

% ENV.Avg_Conv.e_lim_mat = P_diff_mat_in;
% Conv_energy_rating_partition_mat = zeros(CL.Stat.Bat_num + CL.Stat.Bus_num, CL.Stat.Bat_num + CL.Stat.Bus_num);

CL = Update_Bat_stat(Bat_info, CL);
Bat_energy_rating = zeros(CL.Stat.Bat_num,1);
for i = 1:CL.Stat.Bat_num
    Bat_energy_rating(i) = CL.Bat{i}.volt*CL.Bat{i}.qlim;
end

eta = 1; % assuming the converter efficiency is 100% 
Node_num = CL.Stat.Bat_num + CL.Stat.Bus_num + CL.Stat.Cap_num - 1;

P_diff_avg_mat = CL.Conn.diff_2 | CL.Conn.diff_2';    % the interconnection matrix must be symmetric
P_direct_mat = CL.Conn.direct | CL.Conn.direct';           % direct energy dilivery   
P_lsv_mat = CL.Conn.diff_1 | CL.Conn.diff_1';
%Con_direct_num = sum(sum(P_direct_mat))/2;           % count the number of direct connections
%Con_diff_num = sum(sum(P_diff_avg_mat))/2;           % count the number of differential energy delivery connections
% P_diff_var_mat = zeros(11,11);  % the interconnection matrix for variance converter must be symmetric
%lambda = ENV.Constraint.lambda;  % ratio of processed energy/output energy
P_lsv_lim = ENV.Fir_Conv.e_lim_single(k2);  % 0, energy limits of direc-energy-processing converters
P_avg_lim = ENV.Sec_Conv.e_lim_single(k2);
%% Construct the optimization problems and solve the linear programming problems 
P = optimvar('P',Node_num, Node_num);    % dimensional energy transfer variables, primary layer
% P_lsv = optimvar('P_lsv',Node_num, Node_num);   
% Xu = optimvar('Xu',Con_diff_num,'LowerBound',zeros(Con_diff_num,1));  % 3*1 dimensional non negative slack variables
bus_qlim = optimvar('bus_qlim',size(CL.Bus,2),1);     % bus current capability
bus_volt = [];
for i = 1:size(CL.Bus,2)
    CL.Bus{i}.qlim = bus_qlim(i);
    CL.Bus{i}.volt = 0;
    for  j = 1:size(CL.Bat,2)
       CL.Bus{i}.volt = CL.Bus{i}.volt + CL.Bat{j}.volt;
    end
    bus_volt(i) = CL.Bus{i}.volt;
end
prob = optimproblem('Objective',sum(bus_qlim.*bus_volt),'ObjectiveSense','min');  % minimize the summations
% Construct the constraint for the slack energy process ratio 
%energy_ratio_con = optimconstr(1);
% energy_ratio_con = Xu'*ones(Con_diff_num,1) <= -ENV.Constraint.lambda*sum(bus_qlim.*bus_volt); 
%prob.Constraints.energy_ratio_con = energy_ratio_con;

% Construct the constraint for direct energy processing links  
% note that there is not any direct energy processing links 
topo_direct_con = optimconstr(0);
for i = 1:size(CL.Bat,2)
    for j = 1:size(CL.Bus,2)
        if (P_direct_mat(CL.Bat{i}.ind,CL.Bus{j}.ind) == 1)
            topo_direct_con(end+1) = P(CL.Bat{i}.ind,CL.Bus{j}.ind) == - CL.Bat{i}.volt*CL.Bus{j}.qlim;
        end
        if (P_direct_mat(CL.Bus{j}.ind, CL.Bat{i}.ind) == 1)
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

% Construct the constriant for full energy processing links
topo_diff_avg_con = optimconstr(0);
for i = 1:Node_num % 12
    for j = i:Node_num
        if (P_direct_mat(i,j) == 0)
            if (P_diff_avg_mat(i,j) == 1)
                topo_diff_avg_con(end+1) = P(i,j) == -P(j,i);
            % else
            %     topo_diff_avg_con(end+1) = P(i,j) == 0;
            %     topo_diff_avg_con(end+1) = P(j,i) == 0; 
            end
        end
    end
end
 prob.Constraints.topo_diff_avg_con = topo_diff_avg_con;
 
topo_lsv_con = optimconstr(0);
for i = 1:Node_num
    for j = i:Node_num
        if (P_direct_mat(i,j) == 0 && P_diff_avg_mat(i,j) == 0)
            if (P_lsv_mat(i,j) == 1)
                topo_lsv_con(end+1) = P(i,j) == -P(j,i);
            else
                topo_lsv_con(end+1) = P(i,j) == 0;
                topo_lsv_con(end+1) = P(j,i) == 0; 
            end
        end
    end
end
 prob.Constraints.topo_lsv_con = topo_lsv_con;

% % Construct the constriant for differential energy processing links, primary
avg_conv_con = optimconstr(0);
temp_ct = 1;
for i = 1:Node_num
    for j = i:Node_num
        if (P_diff_avg_mat(i,j) == 1)
            avg_conv_con(end+1) = P(i,j) <=  P_avg_lim;
            avg_conv_con(end+1) = -P_avg_lim <= P(i,j);
%            Conv_energy_rating_partition_mat(i,j) = ENV.Avg_Conv.e_lim(temp_ct);
        end
    end
end
prob.Constraints.avg_conv_con = avg_conv_con;


% % (average converters)
% topo_diff_avg_con = optimconstr(0);
% for i = 1:Node_num
%     for j = i:Node_num
%         if (P_fpp_mat(i,j) == 0)
%             if (P_diff_avg_mat(i,j) == 1)
%                 topo_diff_avg_con(end+1) = P(i,j) == -P(j,i);
%             else
%                 topo_diff_avg_con(end+1) = P(i,j) == 0;
%                 topo_diff_avg_con(end+1) = P(j,i) == 0; 
%             end
%         end
%     end
% end
%  prob.Constraints.topo_diff_avg_con = topo_diff_avg_con;
% 
% % Construct the constriant for differential energy processing links,
% % secondary (var converters)
% topo_diff_var_con = optimconstr(0);
% for i = 1:Node_num
%     for j = i:Node_num
%         if (P_diff_var_mat(i,j) == 1)
%             topo_diff_var_con(end+1) = P_var(i,j) == -P_var(j,i);
%         else
%             topo_diff_var_con(end+1) = P_var(i,j) == 0;
%             topo_diff_var_con(end+1) = P_var(j,i) == 0; 
%         end
%     end
% end
%  prob.Constraints.topo_diff_var_con = topo_diff_var_con;

% % Construct the primary/average converter energy limits
% avg_conv_con = optimconstr(0);
% temp_ct = 1;
% for i = 1:Node_num
%     for j = i:Node_num
%         if (P_diff_avg_mat(i,j) == 1)
%             avg_conv_con(end+1) = P(i,j) <=  P_diff_rating_mat_in(i,j);
%             avg_conv_con(end+1) = -P_diff_rating_mat_in(i,j) <= P(i,j);
% %            Conv_energy_rating_partition_mat(i,j) = ENV.Avg_Conv.e_lim(temp_ct);
%         end
%     end
% end
% prob.Constraints.avg_conv_con = avg_conv_con; 

% Construct the lsv-ppp converter energy limits
lsv_conv_con = optimconstr(0);
for i = 1:Node_num
    for j = i:Node_num
        if (P_lsv_mat(i,j) == 1)
            lsv_conv_con(end+1) = P(i,j) <= P_lsv_lim;
            lsv_conv_con(end+1) = -P_lsv_lim <= P(i,j);
%            Conv_energy_rating_partition_mat(i,j) = ENV.Avg_Conv.e_lim(temp_ct);
        end
    end
end
prob.Constraints.lsv_conv_con = lsv_conv_con;

% % Construct the primary/average converter energy limits
% avg_conv_con = optimconstr(0);
% temp_ct = 1;
% for i = 1:Node_num
%     for j = i:Node_num
%         if (P_diff_avg_mat(i,j) == 1)
%             avg_conv_con(end+1) = P(i,j) <=  P_diff_rating_mat_in(i,j);
%             avg_conv_con(end+1) = -P_diff_rating_mat_in(i,j) <= P(i,j);
% %            Conv_energy_rating_partition_mat(i,j) = ENV.Avg_Conv.e_lim(temp_ct);
%         end
%     end
% end
% prob.Constraints.avg_conv_con = avg_conv_con;

% % Construct the energy limit for secondary/variance energy converters
% var_conv_con = optimconstr(0);
% for i = 1:Node_num
%     for j = i:Node_num
%         if (P_diff_var_mat(i,j) == 1)
%             var_conv_con(end+1) = P_var(i,j) <= P_var_lim;
%             var_conv_con(end+1) = -P_var_lim <=  P_var(i,j);
%         end
%     end
% end
% prob.Constraints.var_conv_con = var_conv_con;

% Construct the constriant for energy conservation
energy_bus_con = optimconstr(0);
for i = 1:size(CL.Bus,2)
    energy_bus_con(end+1) = sum(P(CL.Bus{i}.ind,:)) == CL.Bus{i}.volt*CL.Bus{i}.qlim;
end
prob.Constraints.energy_bus_con = energy_bus_con;

energy_cap_con = optimconstr(0);
for i = 1:2
    energy_cap_con(end+1) = sum(P(CL.Cap{i}.ind,:)) == 0;
    % energy_cap_con(end+1) = sum(P_lsv(CL.Cap{i}.ind,:)) == 0;
end
prob.Constraints.energy_cap_con = energy_cap_con;

energy_bat_con = optimconstr(0);
for i = 1:size(CL.Bat,2)
    energy_bat_con(end+1) = CL.Bat{i}.volt*CL.Bat{i}.qlim >= sum(P(CL.Bat{i}.ind,:)) ;
    energy_bat_con(end+1) = sum(P(CL.Bat{i}.ind,:))  >= -CL.Bat{i}.volt*CL.Bat{i}.qlim;
end
prob.Constraints.energy_bat_con = energy_bat_con;   


problem = prob2struct(prob);
[sol,fval,exitflag,output] = linprog(problem);

if(exitflag == 1)
    OptRes.sol = sol;
    OptRes.fval = fval;
    OptRes.P_diff_avg = P_diff_avg_mat;
    % OptRes.P_diff_var = P_diff_var_mat;
    OptRes.P_direct = P_direct_mat;
    OptRes.P_lsv = P_lsv_mat;
    % OptRes.Conv_energy_rating_partition_mat = zeros(11,11);
    
    Total_energy_flow = reshape(OptRes.sol(1:Node_num*Node_num),[Node_num,Node_num]);
    [OptRes.Avg_energy_flow,OptRes.Lsv_energy_flow, OptRes.Direct_power_flow] = func_energy_flow_decomp(Total_energy_flow, CL);
    % OptRes.Var_energy_flow = reshape(OptRes.sol(Node_num*Node_num+1:2*Node_num*Node_num),[Node_num,Node_num]);
    % OptRes.Var_energy_flow = zeros(11,11);
    % OptRes.Var_energy_process = sum(sum(abs(OptRes.Var_energy_flow)))/2;
    OptRes.Avg_energy_process = sum(sum(abs(OptRes.Avg_energy_flow)))/2;
    OptRes.Lsv_energy_process = sum(sum(abs(OptRes.Lsv_energy_flow)))/2;
    OptRes.Total_energy_process = OptRes.Lsv_energy_process + OptRes.Avg_energy_process;% + OptRes.Var_energy_process + 
    
    OptRes.ComponentInfo = CL;

    OptRes.Maximum_output_energy = fval;
    for i = 1: CL.Stat.Bat_num
        OptRes.BatOutputEnergy(i) = sum((OptRes.Lsv_energy_flow(CL.Bat{i}.ind,:))) ...
         + sum((OptRes.Direct_power_flow(CL.Bat{i}.ind,:))) + sum((OptRes.Avg_energy_flow(CL.Bat{i}.ind,:)));
    end
    % OptRes.Var_conv_energy_rating = ones(1,Var_conv_num)*P_var_lim;
    OptRes.Bat_energy_rating = Bat_energy_rating;
    
    % if (sol(end) == sol(2*Node_num*Node_num + 1))
    %     OptRes.Output_q_max = sol(2*Node_num*Node_num + 1:end);
    % else
    %     error('Error in interpreting results in dc energy flow optimization of two lay converters\n');
    % end
else
    CL.Bat{i}.qlim_var 
    ENV.Sec_Conv.e_lim_single
    ENV.Fir_Conv.e_lim_single
    error('Error in dc energy flow optimization on lsv-ppp converters\n');
end

function y = Update_Bat_stat(Bat_Info,CL)
    for l = 1:CL.Stat.Bat_num
        %CL.Bat{l}.qlim = mvnrnd([CL.Bat{l}.qlim_mu,CL.Bat{l}.volt_mu], diag([CL.Bat{l}.qlim_var,CL.Bat{l}.volt_var]),1);
         CL.Bat{l}.qlim = Bat_Info(l);
    end
    CL.Cap{1}.volt = CL.Bus{1}.volt - CL.Bat{1}.volt - CL.Bat{2}.volt - CL.Bat{3}.volt...
    - CL.Bat{4}.volt - CL.Bat{5}.volt - CL.Bat{6}.volt...
    - CL.Bat{7}.volt - CL.Bat{8}.volt - CL.Bat{9}.volt;     % Decide the capactior votlages
    y = CL;
end

function [Avg_energy_flow,Lsv_energy_flow, Direct_energy_flow] = func_energy_flow_decomp(Total_energy_flow,CL)
    Avg_energy_flow = zeros(size(Total_energy_flow,1),size(Total_energy_flow,2));
    Lsv_energy_flow = zeros(size(Total_energy_flow,1),size(Total_energy_flow,2));
    Direct_energy_flow = zeros(size(Total_energy_flow,1),size(Total_energy_flow,2));
    for ic = 1:size(Total_energy_flow,1)
        for jc = 1:size(Total_energy_flow,2)
           if ((ic > CL.Cap{1}.ind)||(jc > CL.Cap{1}.ind)) %12 Avg
               Avg_energy_flow(ic,jc) = Total_energy_flow(ic,jc);
           elseif ((ic == CL.Cap{1}.ind)||(jc == CL.Cap{1}.ind)) % 11 Fir
               Lsv_energy_flow(ic,jc) = Total_energy_flow(ic,jc);
           else
               Direct_energy_flow(ic,jc) = Total_energy_flow(ic,jc);
           end
        end
    end 
end
end
