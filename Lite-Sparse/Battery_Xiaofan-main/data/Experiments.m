%% User Inputs
clear all;
clc;

Bat_info_tab(:,1) = [33.6; 16.8;20.82;10.41;13.02];  % battery power 1st column
Bat_info_tab(:,2) = [11.2; 5.6; 6.94; 3.47; 4.34]; % battery energy 2nd column
Bat_info_tab(:,3) = [29.2; 29.5;20.6;20.6;0.1]/100; % percentage ratio 3rd column
RealBatStat = Real_Battery_Stat(Bat_info_tab);

if (~isfile('ENVEXP.mat'))
%    ENVEXP.RealData.Battery.volt = 12;       % Unit: Volt   
    ENVEXP.RealData.Battery.curlim_mu = 12;  % Unit: A  (?A to ?A)
    ENVEXP.RealData.Battery.qlim_mu = 1;     % Unit: A*h  (?Ah to ?Ah)
    ENVEXP.RealData.Battery.curlim_var = Inf;  % Unit: A
    ENVEXP.RealData.Battery.qlim_var = Inf;    % Unit: A*h
    
    ENVEXP.RealData.Battery.num = 9;   

    %ENVEXP.NomData.Battery.volt = 1;       % we assign nominal battery votlage to be 1   
    ENVEXP.NomData.Battery.curlim_mu = 1;  % we assign nominal battery current to be 1   
    ENVEXP.NomData.Battery.qlim_mu = 1;    % we assign nominal battery capacity to be 1

%   ENVEXP.Base.volt = ENVEXP.RealData.Battery.volt/ENVEXP.NomData.Battery.volt; % Unit: Volt
    ENVEXP.Base.qlim = ENVEXP.RealData.Battery.qlim_mu/ENVEXP.NomData.Battery.qlim_mu;  % Unit: A
    ENVEXP.Base.curlim = ENVEXP.RealData.Battery.curlim_mu/ENVEXP.NomData.Battery.curlim_mu;   % Unit: A*h
 %  ENVEXP.Base.power = ENVEXP.Base.volt*ENVEXP.Base.curlim;
 %  ENVEXP.Base.energy = ENVEXP.Base.volt*ENVEXP.Base.qlim;  % Unit: W*h
    ENVEXP.Base.time = ENVEXP.Base.qlim/ENVEXP.Base.curlim; % Unit: h

    ENVEXP.NomData.Battery.curlim_var = ENVEXP.RealData.Battery.curlim_var/ENVEXP.Base.curlim; 
    ENVEXP.NomData.Battery.qlim_var = ENVEXP.RealData.Battery.qlim_var/ENVEXP.Base.qlim;     
    % global CL ENV;
    % User Input begin
    % please always input Baterry first, then capacitor and at last bus
    % 3 batteries, 1 capacitor, single bus structure

    for i = 1: ENVEXP.RealData.Battery.num
        Input{i} =  struct('sn','bat_lpf','ind',i,'volt',1,'capci',10,'ohmres',0.01,...
           'curlim',ENVEXP.NomData.Battery.curlim_mu,'curlim_var',ENVEXP.NomData.Battery.curlim_var,'curlim_mu',ENVEXP.NomData.Battery.curlim_mu,...
           'curcdf',RealBatStat.norm_value.power_dist.cdf,'curinvcdf', RealBatStat.norm_value.power_dist.invcdf,...
           'res',0.001,'res_mu',0.003,'res_var',0.001/3,...
           'qlim',ENVEXP.NomData.Battery.qlim_mu,'qlim_mu',ENVEXP.NomData.Battery.qlim_mu,'qlim_var',ENVEXP.NomData.Battery.qlim_var,...
           'qcdf',RealBatStat.norm_value.energy_dist.cdf,'qinvcdf', RealBatStat.norm_value.energy_dist.invcdf,...
           'eloss_coeff',0.12,'eloss_coeff_mu',0.12,'eloss_coeff_var',0.012);  
    end
    Input{end+1} = struct('sn','bus1','ind',1,'volt',1,'curlim',-1);
    Input{end+1} = struct('sn','cap1','ind',1,'volt',10,'capci',10000,'ohmres',0.001,'curlim',100);
    % simulation parameters
    %ENV = struct('Constraint',[],'Sweep',[],'Avg_Conv',[],'Var_Conv',[]);   % Initialize the environment
    %ENV.Constraint.lambda = 0.3;
    %ENV.Sweep.Stat.Conv = 5;    % number of var converter power limit to sweep
    %ENV.Sweep.Stat.Bat = 4;     % number of battery types to sweep 
    %ENV.Sweep.Conv.p_lim = linspace(0.001,0.2,ENV.Sweep.Stat.Conv);
    %ENV.Sweep.Conv.e_lim = linspace(0.005,0.05,ENV.Sweep.Stat.Conv);
%     ENVEXP.Avg_Conv.Num = 4;       % the upper num lim of averaging converter
%     %NVEXP.Avg_Conv.ActNum = 3;    % the actural num of averaging converter can be <= the num lim of averaging converter, becuase some of them does not deliver power
%     ENVEXP.Avg_Conv.p_lim_vec = 100*ones(1,ENVEXP.Avg_Conv.Num);  % the upper power lim of averaging converter
%     ENVEXP.Avg_Conv.e_lim_vec = 100*ones(1,ENVEXP.Avg_Conv.Num);  % the upper power lim of averaging converter
%     ENVEXP.Var_Conv.Num = 9;       % the upper num lim of var converter
%     ENVEXP.Var_Conv.p_lim_singlevar = 100; % the upper power lim of var converter, all var converters have the same power limit
%     ENVEXP.Var_Conv.e_lim_singlevar = 100; % the upper power lim of var converter, all var converters have the same power limit
    %ENV.Var_Conv.MC_trial = 2;  % takes 1 hours
    % Categorize the Component lists into Batteries, Cpacitors and Buses
    CL = struct('Bat',[],'Cap',[],'Bus',[],'Conn',[],'Stat',[]);   % Initialize the compoennt list struct
    for i = 1:size(Input,2)
        Input{i}.ind = i;  % updated order is bat, bus and cap
        if strcmp(Input{i}.sn(1:3),'bat')
            CL.Bat{end+1} = Input{i};
        elseif strcmp(Input{i}.sn(1:3),'bus')
            CL.Bus{end+1} = Input{i};
        elseif strcmp(Input{i}.sn(1:3),'cap')
            CL.Cap{end+1} = Input{i};
        else
            fprintf('Error in input\n');
            pause;
        end
    end
    CL.Stat.Bat_num = size(CL.Bat,2);
    CL.Stat.Cap_num = size(CL.Cap,2);
    CL.Stat.Bus_num = size(CL.Bus,2);
    %CL.Stat.Node_num = size(CL.Bat,2) + size(CL.Cap,2) + size(CL.Bus,2);
    for i = 1:CL.Stat.Cap_num   %Initialize the capactior votlages
        CL.Cap{i}.volt = Inf;  
    end
    % Report the Component lists
    fprintf('%d%s\n%d%s\n%d%s\n',size(CL.Bat,2),' Batteries',size(CL.Cap,2),' SuperCapacitors',size(CL.Bus,2),' Buses');
    % for i = 1:CL.Stat.Bat_num
    %     ENV.Sweep.Bat{i}.curlim_mu = linspace(2,2,ENV.Sweep.Stat.Bat);
    %     ENV.Sweep.Bat{i}.curlim_var = linspace(0.1,0.4,ENV.Sweep.Stat.Bat);
    %     ENV.Sweep.Bat{i}.res_mu = linspace(0,0,ENV.Sweep.Stat.Bat);
    %     ENV.Sweep.Bat{i}.res_var = linspace(0,0,ENV.Sweep.Stat.Bat);
    % end

    P_direct_mat_in = zeros(CL.Stat.Bat_num + CL.Stat.Bus_num, CL.Stat.Bat_num + CL.Stat.Bus_num);  % initialize direct power delivery interconnection matrix for P
    P_diff_mat_in = zeros(CL.Stat.Bat_num + CL.Stat.Bus_num, CL.Stat.Bat_num + CL.Stat.Bus_num);    % initialize differential power delivery interconnection matrix for P  
    % Input the direct power delivery connection
    P_direct_mat_in(CL.Bat{1}.ind, CL.Bus{1}.ind) = 1;
    P_direct_mat_in(CL.Bat{2}.ind, CL.Bus{1}.ind) = 1;
    P_direct_mat_in(CL.Bat{3}.ind, CL.Bus{1}.ind) = 1;
    P_direct_mat_in(CL.Bat{4}.ind, CL.Bus{1}.ind) = 1;
    P_direct_mat_in(CL.Bat{5}.ind, CL.Bus{1}.ind) = 1;
    P_direct_mat_in(CL.Bat{6}.ind, CL.Bus{1}.ind) = 1;
    P_direct_mat_in(CL.Bat{7}.ind, CL.Bus{1}.ind) = 1;
    P_direct_mat_in(CL.Bat{8}.ind, CL.Bus{1}.ind) = 1;
    P_direct_mat_in(CL.Bat{9}.ind, CL.Bus{1}.ind) = 1;
    % P_direct_mat_in(CL.Cap{1}.ind, CL.Bus{1}.ind) = 1;
    % Input the direct power delivery connection
    P_diff_var_mat_in = zeros(CL.Stat.Bat_num + CL.Stat.Bus_num, CL.Stat.Bat_num + CL.Stat.Bus_num);
    P_diff_var_mat_in(CL.Bat{1}.ind, CL.Bat{2}.ind) = 1;
    P_diff_var_mat_in(CL.Bat{2}.ind, CL.Bat{3}.ind) = 1;
    P_diff_var_mat_in(CL.Bat{3}.ind, CL.Bat{4}.ind) = 1;
    P_diff_var_mat_in(CL.Bat{4}.ind, CL.Bat{5}.ind) = 1;
    P_diff_var_mat_in(CL.Bat{5}.ind, CL.Bat{6}.ind) = 1;
    P_diff_var_mat_in(CL.Bat{6}.ind, CL.Bat{7}.ind) = 1;
    P_diff_var_mat_in(CL.Bat{7}.ind, CL.Bat{8}.ind) = 1;
    P_diff_var_mat_in(CL.Bat{8}.ind, CL.Bat{9}.ind) = 1;
    % P_diff_var_mat_in(CL.Bat{9}.ind, CL.Cap{1}.ind) = 1;
    % the traditional DPP architectures
    %User Input end

    % Quick Check
    InputCheck(CL,ENVEXP);
    CL.Conn.direct = P_direct_mat_in;
    CL.Conn.diff_var = P_diff_var_mat_in;
    CL.Conn.diff_avg = zeros(CL.Stat.Bat_num + CL.Stat.Bus_num, CL.Stat.Bat_num + CL.Stat.Bus_num);
    CL.Stat.VarConn_num = sum(sum(P_diff_var_mat_in));

% Generate a specific battery statistics
%   ENVEXP.NomData.Battery.curlim = normrnd(ENVEXP.NomData.Battery.curlim_mu, ENVEXP.NomData.Battery.curlim_var, ENVEXP.RealData.Battery.num,1);
    ENVEXP.NomData.Battery.curlim = RealBatStat.norm_value.power_dist.invcdf(rand(ENVEXP.RealData.Battery.num,1));   
    ENVEXP.NomData.Battery.curlim = sort(ENVEXP.NomData.Battery.curlim);  
%   ENVEXP.NomData.Battery.qlim = normrnd(ENVEXP.NomData.Battery.qlim_mu, ENVEXP.NomData.Battery.qlim_var,ENVEXP.RealData.Battery.num,1);
    ENVEXP.NomData.Battery.qlim = RealBatStat.norm_value.energy_dist.invcdf(rand(ENVEXP.RealData.Battery.num,1));  
    ENVEXP.NomData.Battery.qlim = sort(ENVEXP.NomData.Battery.qlim);

% Assign these battery statistitics to CL     
    for i = 1:CL.Stat.Bat_num
        CL.Bat{i}.curlim = ENVEXP.NomData.Battery.curlim(i);
        CL.Bat{i}.qlim = ENVEXP.NomData.Battery.qlim(i);
        CL.Bat_Info.Current(i) = CL.Bat{i}.curlim;
 %      CL.Bat_Info.Power(i) = CL.Bat{i}.curlim*CL.Bat{i}.volt;
        CL.Bat_Info.Charge(i) = CL.Bat{i}.qlim;
 %      CL.Bat_Info.Energy(i) = CL.Bat{i}.qlim*CL.Bat{i}.volt;
    end
    
    ENVEXP.RealData.Battery.curlim = ENVEXP.NomData.Battery.curlim*ENVEXP.Base.curlim;
    ENVEXP.RealData.Battery.qlim = ENVEXP.NomData.Battery.qlim*ENVEXP.Base.qlim;
    
 %   ENVEXP.RealData.Battery.plim = ENVEXP.RealData.Battery.volt*ENVEXP.RealData.Battery.curlim;
 %   ENVEXP.RealData.Battery.elim = ENVEXP.RealData.Battery.volt*ENVEXP.RealData.Battery.qlim;
    
    save('ENVEXP','ENVEXP');
    save('CL','CL');
else

end
load('ENVEXP');
ENVEXP.Avg_Conv.Num = 3;       % the upper num lim of averaging converter
%NVEXP.Avg_Conv.ActNum = 3;    % the actural num of averaging converter can be <= the num lim of averaging converter, becuase some of them does not deliver power
ENVEXP.Avg_Conv.p_lim_vec = 100*ones(1,ENVEXP.Avg_Conv.Num);  % the upper power lim of averaging converter
ENVEXP.Avg_Conv.e_lim_vec = 100*ones(1,ENVEXP.Avg_Conv.Num);  % the upper power lim of averaging converter
ENVEXP.Var_Conv.Num = 8;       % the  num lim of var converter
ENVEXP.Var_Conv.p_lim_singlevar = 100; % the upper power lim of var converter, all var converters have the same power limit
ENVEXP.Var_Conv.e_lim_singlevar = 100; % the upper power lim of var converter, all var converters have the same power limit
% update the votlages vectors
ENVEXP.RealData.Battery.volt = [14, 14, 14, 14, 14, 14, 14, 14, 14];
ENVEXP.Base.volt = 12;       % we assign nominal battery votlage to be 12 
ENVEXP.NomData.Battery.volt = ENVEXP.RealData.Battery.volt./ENVEXP.Base.volt;
ENVEXP.RealData.Battery.plim = ENVEXP.RealData.Battery.volt'.*ENVEXP.RealData.Battery.curlim;
ENVEXP.RealData.Battery.elim = ENVEXP.RealData.Battery.volt'.*ENVEXP.RealData.Battery.qlim;
save('ENVEXP');
%% Task A: Design topology and sizing A by output power optimization
clear;
load('CL');
load('ENVEXP');
ENVA = ENVEXP;
ENVA.Avg_Conv.trial_num = 1e4; % the numer of random searches to find the optimal placement of averaging converter
ENVA.Avg_Conv.partition = 1;   % partition the averaging converter into groups
for i = 1:CL.Stat.Bat_num
    CL.Bat{i}.volt = ENVEXP.NomData.Battery.volt(i);
    CL.Bat_Info.Power(i) = CL.Bat{i}.curlim*CL.Bat{i}.volt;
    CL.Bat_Info.Energy(i) = CL.Bat{i}.qlim*CL.Bat{i}.volt;
end
OptAvg_power = func_average_converter_power_design(CL,ENVA);
% postprocessing: sort the cnverter rating in ascending orders
[OptAvg_power.Conv_power_rating_partition, OptAvg_power.Conv_power_rating_partition_mat] ...
    = func_rating_partition(abs(OptAvg_power.Avg_power_flow), ENVA.Avg_Conv.partition); 
% postprocessing: update the primary converter design results to ENV4
ENVA.Avg_Conv.p_lim_vec = OptAvg_power.Conv_power_rating_partition;
ENVA.Avg_Conv.p_lim_mat = OptAvg_power.Conv_power_rating_partition_mat;
ENVA.Avg_Conv.Num = sum(sum(OptAvg_power.Conv_power_rating_partition_mat > 0))/2;

fprintf('The Battery Power Utilization Ratio is\n');    % ~75%
ENVA.Avg_Conv.bat_power_uratio = abs(OptAvg_power.output_power)/sum(OptAvg_power.Bat_power_rating);
ENVA.Avg_Conv.bat_power_uratio
fprintf('The Power Utilization ratio of the Average Converters\n');
ENVA.Avg_Conv.power_uratio = sum(sum(abs(OptAvg_power.Avg_power_flow)))/sum(sum(ENVA.Avg_Conv.p_lim_mat));
ENVA.Avg_Conv.power_uratio

PreDesignA = OptAvg_power;
save('PreDesignA','PreDesignA');
%save('OptAvg_3p','OptAvg_power');
% ENVA.Sweep.Stat.Bat = 1; 
% for i = 1:CL.Stat.Bat_num
%     ENVA.Sweep.Bat{i}.curlim_mu = ENVA.NomData.Battery.curlim_mu ;
%     ENVA.Sweep.Bat{i}.curlim_var = ENVA.NomData.Battery.curlim_var;
% end
% ENVA.Sweep.Stat.Conv = 5;    % number of var converter power limit to sweep
% ENVA.Sweep.Conv.p_lim = linspace(0.1,1,ENVA.Sweep.Stat.Conv)*sum(ENVA.Avg_Conv.p_lim)/ENVA.Var_Conv.Num;
% ENVA.Var_Conv.MC_trial = 2;  % takes 1 hours
% tic
% DesignA1 = func_var_converter_power_design(CL,ENVA);
% save('DesignA1','DesignA1');
% toc
% func_plot_task0(CL,ENVA,DesignA1);
% fprintf('pick the var converters by lambda please\n')
%pause;
%ENVA.Var_Conv.p_lim = 0.4*sum(ENVA.Avg_Conv.p_lim)/CL.Stat.VarConn_num; 
ENVA.Var_Conv.p_lim_singlevar = 0.2*min(ENVA.Avg_Conv.p_lim_vec);  % select the var converters so that it is 20% of the primary converters
% We purchase the batteries
ENVA.RealData.Battery.volt = [15, 15, 15, 13, 13, 13, 13, 13, 13];
ENVA.Base.volt = 12;       % we assign nominal battery votlage to be 12 
ENVA.NomData.Battery.volt = ENVA.RealData.Battery.volt./ENVA.Base.volt;
for i = 1:CL.Stat.Bat_num
    CL.Bat{i}.volt = ENVA.NomData.Battery.volt(i);
    CL.Bat_Info.Power(i) = CL.Bat{i}.curlim*CL.Bat{i}.volt;
    CL.Bat_Info.Energy(i) = CL.Bat{i}.qlim*CL.Bat{i}.volt;
end
ENVA.RealData.Battery.plim = ENVA.RealData.Battery.volt'.*ENVA.RealData.Battery.curlim;
ENVA.RealData.Battery.elim = ENVA.RealData.Battery.volt'.*ENVA.RealData.Battery.qlim;
% After purchasing the converters and batteries, we reorganize the topology
DesignA = func_two_converter_power_design(CL,ENVA);
% if(size(DesignA_collection,2) > 1)
%     pause;
%     fprintf('Multiple Optimal Points \n');
%     DesignA = DesignA_collection{2};
% else
%     DesignA = DesignA_collection{1};
% end
ENVA.Avg_Conv.p_lim_mat = DesignA.P_diff_avg_mat_rating;
save('DesignA','DesignA');
save('CL');
save('ENVA','ENVA');
%% Task A1 - Use Circuit A for power trajectory (Trajectory 1), design converter flow by power optimization (Exp 1)
clear;
load('CL');
load('ENVA');
load('DesignA','DesignA');

ENVA.RealData.Bus.power_trajectory = -[0 1000 0];   % Unit: Watts
ENVA.RealData.Bus.volt = 130;     % Unit: Volt
ENVA.NomData.Bus.volt = ENVA.RealData.Bus.volt/ENVA.Base.volt; 
%ENVA.RealData.Bus.current_trajectory = ENVA.RealData.Bus.power_trajectory/ENVA.RealData.Bus.volt;   % Unit: Watts
ENVA.NomData.Bus.power_trajectory = ENVA.RealData.Bus.power_trajectory/ENVA.Base.power;  

% Simulations
for i = 1:length(ENVA.NomData.Bus.power_trajectory)
    if (abs(ENVA.NomData.Bus.power_trajectory(i)) >= abs(DesignA.Maximum_output_power))
        CurDist_A1{i}.BESSDiscountFactor = 1;
        CurDist_A1{i}.Nomdata.Bus_Power = -DesignA.Maximum_output_power;
    else
        CurDist_A1{i}.BESSDiscountFactor = abs(ENVA.NomData.Bus.power_trajectory(i)/sum(DesignA.Maximum_output_power));
        CurDist_A1{i}.Nomdata.Bus_Power = ENVA.NomData.Bus.power_trajectory(i);
    end
    CurDist_A1{i}.Nomdata.Avg_Power = sum(DesignA.Avg_power_flow(1:ENVA.RealData.Battery.num,1:ENVA.RealData.Battery.num),2)*CurDist_A1{i}.BESSDiscountFactor;
    CurDist_A1{i}.Nomdata.Var_Power = sum(DesignA.Var_power_flow(1:ENVA.RealData.Battery.num,1:ENVA.RealData.Battery.num),2)*CurDist_A1{i}.BESSDiscountFactor;
    CurDist_A1{i}.Nomdata.Sum_Power = CurDist_A1{i}.Nomdata.Avg_Power + CurDist_A1{i}.Nomdata.Var_Power;
    for j = 1:ENVA.RealData.Battery.num 
        CurDist_A1{i}.Nomdata.Avg_Current(j) = CurDist_A1{i}.Nomdata.Avg_Power(j)./CL.Bat{j}.volt;
        CurDist_A1{i}.Nomdata.Var_Current(j) = CurDist_A1{i}.Nomdata.Var_Power(j)./CL.Bat{j}.volt;
        CurDist_A1{i}.Nomdata.Sum_Current(j) = CurDist_A1{i}.Nomdata.Sum_Power(j)./CL.Bat{j}.volt;
    end
    CurDist_A1{i}.Nomdata.Bus_Current = CurDist_A1{i}.Nomdata.Bus_Power./ENVA.NomData.Bus.volt;
end
for i = 1:length(ENVA.NomData.Bus.power_trajectory)
    CurDist_A1{i}.Realdata.Avg_Current = CurDist_A1{i}.Nomdata.Avg_Current*ENVA.Base.curlim;
    CurDist_A1{i}.Realdata.Var_Current = CurDist_A1{i}.Nomdata.Var_Current*ENVA.Base.curlim;
    CurDist_A1{i}.Realdata.Sum_Current = CurDist_A1{i}.Nomdata.Sum_Current*ENVA.Base.curlim;
    CurDist_A1{i}.Realdata.Bus_Current = CurDist_A1{i}.Nomdata.Bus_Current*ENVA.Base.curlim;
    CurDist_A1{i}.Alpha.Avg = CurDist_A1{i}.Nomdata.Avg_Current/CurDist_A1{i}.Nomdata.Bus_Current;
    CurDist_A1{i}.Alpha.Var = CurDist_A1{i}.Nomdata.Var_Current/CurDist_A1{i}.Nomdata.Bus_Current;
    CurDist_A1{i}.Alpha.Sum = CurDist_A1{i}.Nomdata.Sum_Current/CurDist_A1{i}.Nomdata.Bus_Current;
end
save('CurDist_A1','CurDist_A1');
% Simu_A1 = {};
% for i = 1:length(ENVA.NomData.Bus.power_trajectory)
%     if (abs(ENVA.NomData.Bus.power_trajectory(i)) >= abs(DesignA.Maximum_output_power))
%         Simu_A1{i}.BESSDiscountFactor = 1;
%         Simu_A1{i}.Nomdata.Bus_Power = -DesignA.Maximum_output_power;
%     else
%         Simu_A1{i}.BESSDiscountFactor = abs(ENVA.NomData.Bus.power_trajectory(i)/sum(DesignA.Maximum_output_power));
%         Simu_A1{i}.Nomdata.Bus_Power = ENVA.NomData.Bus.power_trajectory(i);
%     end
%     Simu_A1{i}.Nomdata.Avg_Power = sum(DesignA.Avg_power_flow(1:ENVA.RealData.Battery.num,1:ENVA.RealData.Battery.num),2)*Simu_A1{i}.BESSDiscountFactor;
%     Simu_A1{i}.Nomdata.Var_Power = sum(DesignA.Var_power_flow(1:ENVA.RealData.Battery.num,1:ENVA.RealData.Battery.num),2)*Simu_A1{i}.BESSDiscountFactor;
%     Simu_A1{i}.Nomdata.Sum_Power = Simu_A1{i}.Nomdata.Avg_Power + Simu_A1{i}.Nomdata.Var_Power;
% end
% for i = 1:length(ENVA.NomData.Bus.power_trajectory)
%     Simu_A1{i}.Realdata.Avg_Power = Simu_A1{i}.Nomdata.Avg_Power*ENVA.Base.power;
%     Simu_A1{i}.Realdata.Var_Power = Simu_A1{i}.Nomdata.Var_Power*ENVA.Base.power;
%     Simu_A1{i}.Realdata.Sum_Power = Simu_A1{i}.Nomdata.Sum_Power*ENVA.Base.power;
%     Simu_A1{i}.Realdata.Bus_Power = Simu_A1{i}.Nomdata.Bus_Power*ENVA.Base.power;
% end
% save('Simu_A1','Simu_A1');

% Hardware experiments
for i = 1:length(ENVA.NomData.Bus.power_trajectory)
    if (abs(ENVA.NomData.Bus.power_trajectory(i)) >= abs(DesignA.Maximum_output_power))
        Exp_A1{i}.BESSDiscountFactor = 1;
        Exp_A1{i}.Nomdata.Bus_Power = -DesignA.Maximum_output_power;
    else
        Exp_A1{i}.BESSDiscountFactor = abs(ENVA.NomData.Bus.power_trajectory(i)/sum(DesignA.Maximum_output_power));
        Exp_A1{i}.Nomdata.Bus_Power = ENVA.NomData.Bus.power_trajectory(i);
    end
    Exp_A1{i}.Nomdata.Avg_Power = DesignA.Avg_power_flow(1:ENVA.RealData.Battery.num,1:ENVA.RealData.Battery.num)*Exp_A1{i}.BESSDiscountFactor;
    Exp_A1{i}.Nomdata.Var_Power = DesignA.Var_power_flow(1:ENVA.RealData.Battery.num,1:ENVA.RealData.Battery.num)*Exp_A1{i}.BESSDiscountFactor;
end   

for i = 1:length(ENVA.NomData.Bus.power_trajectory)
    Exp_A1{i}.Realdata.Avg_Power = Exp_A1{i}.Nomdata.Avg_Power*ENVA.Base.power;
    Exp_A1{i}.Realdata.Var_Power = Exp_A1{i}.Nomdata.Var_Power*ENVA.Base.power;
    Exp_A1{i}.Realdata.Bus_Power = Exp_A1{i}.Nomdata.Bus_Power*ENVA.Base.power;
end
save('Exp_A1','Exp_A1');

save('ENVA','ENVA');

%% Task A2 - Use Circuit A for energy trajectory (Trajectory 2), design converter flow by energy optimization (Exp 2)
clear;
load('CL');
load('ENVA');
load('DesignA');

ENVA.RealData.Bus.energy_trajectory = - [180 180 180];  % Unit: Watts, corresponds to 1.5A output current 
ENVA.RealData.Time.discharge = 3; % Unit: h
ENVA.NomData.Time.discharge = ENVA.RealData.Time.discharge/ENVA.Base.time;

% (Resolved) importantly note that the topology should come from the power optimization
Temp.Design1 = func_dc_energyflow_avg_layer(ENVA.Avg_Conv.p_lim_mat*100, CL, ENVA);
ENVA.Avg_Conv.e_lim_mat = abs(Temp.Design1.Avg_energy_flow);
ENVA.Avg_Conv.e_lim_vec = nonzeros(triu(ENVA.Avg_Conv.e_lim_mat));
ENVA.NomData.Bus.energy_trajectory = ENVA.RealData.Bus.energy_trajectory/ENVA.Base.power;  

fprintf('The battery utilization is\n');    % ~75%
ENVA.Avg_Conv.bat_energy_uratio = abs(Temp.Design1.bus_elim)/sum(CL.Bat_Info.Energy);
ENVA.Avg_Conv.bat_energy_uratio
fprintf('The Energy Utilization ratio of the Average Converters\n');  % ~100%
ENVA.Avg_Conv.energy_uratio = sum(sum(abs(Temp.Design1.Avg_energy_flow)))/sum(sum(ENVA.Avg_Conv.e_lim_mat));
ENVA.Avg_Conv.energy_uratio
%ENVA.Var_Conv.e_lim = 0.4*sum(ENVA.Avg_Conv.e_lim)/CL.Stat.VarConn_num;  % select the var converters using lambda = 0.5
ENVA.Var_Conv.e_lim_singlevar = 0.2*min(ENVA.Avg_Conv.e_lim_vec);  % select the var converters so that it is 20% of the primary converters
DesignA2 = func_dc_energyflow_two_layers(ENVA.Avg_Conv.e_lim_mat, CL.Bat_Info.Charge, CL, ENVA);
save('DesignA2', 'DesignA2');
for i = 1:length(ENVA.NomData.Bus.energy_trajectory)
    if (abs(ENVA.NomData.Bus.energy_trajectory(i)) >= abs(DesignA.Maximum_output_power))
        %CurDist_A2{i}.BESSDiscountFactor = 1;
        CurDist_A2{i}.Nomdata.Bus_Power = - DesignA.Maximum_output_power;
        error('Battery Overdischarge');
    else
        %CurDist_A2{i}.BESSDiscountFactor = abs(ENVA.NomData.Bus.energy_trajectory(i)/sum(DesignA.Maximum_output_power));
        CurDist_A2{i}.Nomdata.Bus_power = ENVA.NomData.Bus.energy_trajectory(i);
        
 %       CurDist_A2{i}.Nomdata.FB_time = DesignA2.Maximum_output_energy/ENVA.NomData.Bus.energy_trajectory(i);
    end
    CurDist_A2{i}.Nomdata.Avg_power = sum(DesignA2.Avg_energy_flow(1:ENVA.RealData.Battery.num,1:ENVA.RealData.Battery.num),2)...
                                        /DesignA2.Maximum_output_energy...
                                        *(ENVA.NomData.Bus.energy_trajectory(i));
    CurDist_A2{i}.Nomdata.Var_power = sum(DesignA2.Var_energy_flow(1:ENVA.RealData.Battery.num,1:ENVA.RealData.Battery.num),2)...
                                        /DesignA2.Maximum_output_energy...
                                        *(ENVA.NomData.Bus.energy_trajectory(i));
    CurDist_A2{i}.Nomdata.Sum_power = CurDist_A2{i}.Nomdata.Avg_power + CurDist_A2{i}.Nomdata.Var_power;
    for j = 1:ENVA.RealData.Battery.num 
        CurDist_A2{i}.Nomdata.Avg_Current(j) = CurDist_A2{i}.Nomdata.Avg_power(j)./CL.Bat{j}.volt;
        CurDist_A2{i}.Nomdata.Var_Current(j) = CurDist_A2{i}.Nomdata.Var_power(j)./CL.Bat{j}.volt;
        CurDist_A2{i}.Nomdata.Sum_Current(j) = CurDist_A2{i}.Nomdata.Sum_power(j)./CL.Bat{j}.volt;
    end
    CurDist_A2{i}.Nomdata.Bus_Current = CurDist_A2{i}.Nomdata.Bus_power./ENVA.NomData.Bus.volt;
end
for i = 1:length(ENVA.NomData.Bus.energy_trajectory)
    CurDist_A2{i}.Realdata.Avg_Current = CurDist_A2{i}.Nomdata.Avg_Current*ENVA.Base.curlim;
    CurDist_A2{i}.Realdata.Var_Current = CurDist_A2{i}.Nomdata.Var_Current*ENVA.Base.curlim;
    CurDist_A2{i}.Realdata.Sum_Current = CurDist_A2{i}.Nomdata.Sum_Current*ENVA.Base.curlim;
    CurDist_A2{i}.Realdata.Bus_Current = CurDist_A2{i}.Nomdata.Bus_Current*ENVA.Base.curlim;
    CurDist_A2{i}.Alpha.Avg = CurDist_A2{i}.Nomdata.Avg_Current/CurDist_A2{i}.Nomdata.Bus_Current;
    CurDist_A2{i}.Alpha.Var = CurDist_A2{i}.Nomdata.Var_Current/CurDist_A2{i}.Nomdata.Bus_Current;
    CurDist_A2{i}.Alpha.Sum = CurDist_A2{i}.Nomdata.Sum_Current/CurDist_A2{i}.Nomdata.Bus_Current;
end
save('CurDist_A2','CurDist_A2');

save('ENVA','ENVA');

%% Task B: design topology and sizing B by energy optimization 
clear;
load('ENVEXP');
load('CL');
ENVB = ENVEXP;
ENVB.Avg_Conv.trial_num = 1e4; % the numer of random searches to find the optimal placement of averaging converter
ENVB.Avg_Conv.partition = 1;   % partition the averaging converter into groups
for i = 1:CL.Stat.Bat_num
    CL.Bat{i}.volt = ENVB.NomData.Battery.volt(i);
    CL.Bat_Info.Power(i) = CL.Bat{i}.curlim*CL.Bat{i}.volt;
    CL.Bat_Info.Energy(i) = CL.Bat{i}.qlim*CL.Bat{i}.volt;
end
OptAvg_energy = func_average_converter_energy_design(CL,ENVB);
% postprocessing: sort the cnverter rating in ascending orders
[OptAvg_energy.Conv_energy_rating_partition, OptAvg_energy.Conv_energy_rating_partition_mat] ...
    = func_rating_partition(abs(OptAvg_energy.Avg_energy_flow), ENVB.Avg_Conv.partition); 
% postprocessing: update the primary converter design results to ENV4
ENVB.Avg_Conv.e_lim_vec = OptAvg_energy.Conv_energy_rating_partition;
ENVB.Avg_Conv.e_lim_mat = OptAvg_energy.Conv_energy_rating_partition_mat;
ENVB.Avg_Conv.Num = sum(sum(OptAvg_energy.Conv_energy_rating_partition_mat > 0))/2;

fprintf('The Battery energy Utilization Ratio is\n');  
ENVB.Avg_Conv.bat_energy_uratio = abs(OptAvg_energy.output_energy)/sum(OptAvg_energy.Bat_energy_rating);
ENVB.Avg_Conv.bat_energy_uratio
fprintf('The energy Utilization ratio of the Average Converters\n');
ENVB.Avg_Conv.energy_uratio = sum(sum(abs(OptAvg_energy.Avg_energy_flow)))/sum(sum(ENVB.Avg_Conv.e_lim_mat));
ENVB.Avg_Conv.energy_uratio
PreDesignB = OptAvg_energy;
save('PreDesignB','PreDesignB');
% Design the varaiance converters
%save('OptAvg_3p','OptAvg_energy');
% ENVB.Sweep.Stat.Bat = 1; 
% for i = 1:CL.Stat.Bat_num
%     ENVB.Sweep.Bat{i}.qlim_mu = ENVB.NomData.Battery.qlim_mu ;
%     ENVB.Sweep.Bat{i}.qlim_var = ENVB.NomData.Battery.qlim_var;
% end
% ENVB.Sweep.Stat.Conv = 5;    % number of var converter energy limit to sweep
% ENVB.Sweep.Conv.e_lim = linspace(0.1,1,ENVB.Sweep.Stat.Conv)*sum(ENVB.Avg_Conv.e_lim)/ENVB.Var_Conv.Num;
% ENVB.Var_Conv.MC_trial = 2;  % takes 1 hours
% tic
% DesignB1 = func_var_converter_energy_design(CL,ENVB);
% save('DesignB1','DesignB1');
% toc
% func_plot_task4(CL,ENVB,DesignB1);
% fprintf('pick the var converters by lambda\n')
% pause;
% We purchase the batteries
ENVB.RealData.Battery.volt = [15, 15, 15, 13, 13, 13, 13, 13, 13];
ENVB.Base.volt = 12;       % we assign nominal battery votlage to be 12 
ENVB.NomData.Battery.volt = ENVB.RealData.Battery.volt./ENVB.Base.volt;
for i = 1:CL.Stat.Bat_num
    CL.Bat{i}.volt = ENVB.NomData.Battery.volt(i);
    CL.Bat_Info.Power(i) = CL.Bat{i}.curlim*CL.Bat{i}.volt;
    CL.Bat_Info.Energy(i) = CL.Bat{i}.qlim*CL.Bat{i}.volt;
end
ENVB.RealData.Battery.plim = ENVB.RealData.Battery.volt'.*ENVB.RealData.Battery.curlim;
ENVB.RealData.Battery.elim = ENVB.RealData.Battery.volt'.*ENVB.RealData.Battery.qlim;

ENVB.Var_Conv.e_lim_singlevar = 0.2*min(ENVB.Avg_Conv.e_lim_vec);  % select the var converters to be 20% of the primary converters
DesignB = func_two_converter_energy_design(CL,ENVB);
% if(size(DesignB_collection,2) > 1)
%     pause;
%     pfrintf('Multiple Optimal Points \n')
%     DesignB = DesignB_collection{2};
% else
%     DesignB = DesignB_collection{1};
% end
ENVB.Avg_Conv.e_lim_mat = DesignB.P_diff_avg_mat_rating;
save('DesignB','DesignB');
% Power Converter Sizing for 3 hours discharging
ENVB.Avg_Conv.p_lim_vec = ENVB.Avg_Conv.e_lim_vec/mean(ENVEXP.NomData.Battery.volt)/ENVEXP.NomData.Battery.qlim_mu...
                       *mean(ENVEXP.NomData.Battery.volt)*ENVEXP.NomData.Battery.curlim_mu/ENVEXP.NomData.Time.discharge;
                   
ENVB.Avg_Conv.p_lim_mat = ENVB.Avg_Conv.e_lim_mat/mean(ENVEXP.NomData.Battery.volt)/ENVEXP.NomData.Battery.qlim_mu...
                       *mean(ENVEXP.NomData.Battery.volt)*ENVEXP.NomData.Battery.curlim_mu/ENVEXP.NomData.Time.discharge;
                   
ENVB.Var_Conv.p_lim_singlevar = ENVB.Var_Conv.e_lim_singlevar/mean(ENVEXP.NomData.Battery.volt)/ENVEXP.NomData.Battery.qlim_mu...
                       *mean(ENVEXP.NomData.Battery.volt)*ENVEXP.NomData.Battery.curlim_mu/ENVEXP.NomData.Time.discharge;
save('ENVB','ENVB');
%% Task B1 - Use Circuit B for power trajectory (Trajectory 1), design converter flow by power optimization  (Exp 3)
clear;
load('CL');
load('ENVB');
DesignB1 = func_dc_powerflow_two_layers(ENVB.Avg_Conv.p_lim_mat, CL.Bat_Info.Current, CL, ENVB);
save('DesignB1','DesignB1');
for i = 1:length(ENVB.NomData.Bus.power_trajectory)
    if (abs(ENVB.NomData.Bus.power_trajectory(i)) >= abs(DesignB1.Maximum_output_power))
        CurDist_B1{i}.BESSDiscountFactor = 1;
        CurDist_B1{i}.Nomdata.Bus_Power = -DesignB1.Maximum_output_power;
    else
        CurDist_B1{i}.BESSDiscountFactor = abs(ENVB.NomData.Bus.power_trajectory(i)/abs(DesignB1.Maximum_output_power));
        CurDist_B1{i}.Nomdata.Bus_Power = ENVB.NomData.Bus.power_trajectory(i);
    end
    CurDist_B1{i}.Nomdata.Avg_Power = sum(DesignB1.Avg_power_flow(1:ENVB.RealData.Battery.num,1:ENVB.RealData.Battery.num),2)*CurDist_B1{i}.BESSDiscountFactor;
    CurDist_B1{i}.Nomdata.Var_Power = sum(DesignB1.Var_power_flow(1:ENVB.RealData.Battery.num,1:ENVB.RealData.Battery.num),2)*CurDist_B1{i}.BESSDiscountFactor;
    CurDist_B1{i}.Nomdata.Sum_Power = CurDist_B1{i}.Nomdata.Avg_Power + CurDist_B1{i}.Nomdata.Var_Power;
    for j = 1:ENVB.RealData.Battery.num 
        CurDist_B1{i}.Nomdata.Avg_Current(j) = CurDist_B1{i}.Nomdata.Avg_Power(j)./CL.Bat{j}.volt;
        CurDist_B1{i}.Nomdata.Var_Current(j) = CurDist_B1{i}.Nomdata.Var_Power(j)./CL.Bat{j}.volt;
        CurDist_B1{i}.Nomdata.Sum_Current(j) = CurDist_B1{i}.Nomdata.Sum_Power(j)./CL.Bat{j}.volt;
    end
    CurDist_B1{i}.Nomdata.Bus_Current = CurDist_B1{i}.Nomdata.Bus_Power./ENVB.NomData.Bus.volt;
end
for i = 1:length(ENVB.NomData.Bus.power_trajectory)
    CurDist_B1{i}.Realdata.Avg_Current = CurDist_B1{i}.Nomdata.Avg_Current*ENVB.Base.curlim;
    CurDist_B1{i}.Realdata.Var_Current = CurDist_B1{i}.Nomdata.Var_Current*ENVB.Base.curlim;
    CurDist_B1{i}.Realdata.Sum_Current = CurDist_B1{i}.Nomdata.Sum_Current*ENVB.Base.curlim;
    CurDist_B1{i}.Realdata.Bus_Current = CurDist_B1{i}.Nomdata.Bus_Current*ENVB.Base.curlim;
    CurDist_B1{i}.Alpha.Avg = CurDist_B1{i}.Nomdata.Avg_Current/CurDist_B1{i}.Nomdata.Bus_Current;
    CurDist_B1{i}.Alpha.Var = CurDist_B1{i}.Nomdata.Var_Current/CurDist_B1{i}.Nomdata.Bus_Current;
    CurDist_B1{i}.Alpha.Sum = CurDist_B1{i}.Nomdata.Sum_Current/CurDist_B1{i}.Nomdata.Bus_Current;
end
save('CurDist_B1','CurDist_B1');
%% Task B2 - Use Circuit B for energy trajectory (Trajectory 2), design converter flow by energy optimization (Exp 4)
clear;
load('CL');
load('ENVB');
load('DesignB');   % Use the energy design in the previous sections
load('DesignB1');  % need to take use of power limitation of BESS in Task B1
%DesignB2 = func_dc_energyflow_two_layers(ENVB.Avg_Conv.e_lim_mat, CL.Bat_Info.Charge, CL, ENVB);
% save('DesignB2', 'DesignB2');
% save('ENVB','ENVB');
ENVB.RealData.Bus.energy_trajectory = -[0 30 0];   % Unit: Watts
ENVB.RealData.Bus.volt = 130;     % Unit: Volt
ENVB.NomData.Bus.volt = ENVB.RealData.Bus.volt/ENVB.Base.volt; 
%ENVB.RealData.Bus.current_trajectory = ENVB.RealData.Bus.power_trajectory/ENVB.RealData.Bus.volt;   % Unit: Watts
ENVB.NomData.Bus.energy_trajectory = ENVB.RealData.Bus.energy_trajectory/ENVB.Base.power;  

for i = 1:length(ENVB.NomData.Bus.energy_trajectory)
    if (abs(ENVB.NomData.Bus.energy_trajectory(i)) >= abs(DesignB1.Maximum_output_power))
        %CurDist_B2{i}.BESSDiscountFactor = 1;
        CurDist_B2{i}.Nomdata.Bus_Power = DesignB1.Maximum_output_power;
        error('BESS Over Discharged');
    else
        %CurDist_B2{i}.BESSDiscountFactor = abs(ENVB.NomData.Bus.energy_trajectory(i)/DesignB1.Maximum_output_power);
    
         %CurDist_B2{i}.Nomdata.Bus_power = -DesignB.Maximum_output_energy/ENVB.NomData.Time.discharge;
        CurDist_B2{i}.Nomdata.Bus_power = ENVB.NomData.Bus.energy_trajectory(i);
        
 %       CurDist_B2{i}.Nomdata.FB_time = DesignB2.Maximum_output_energy/ENVB.NomData.Bus.energy_trajectory(i);
    end
    CurDist_B2{i}.Nomdata.Avg_power = sum(DesignB.Avg_energy_flow(1:ENVB.RealData.Battery.num,1:ENVB.RealData.Battery.num),2)...
                                        /(-DesignB.Maximum_output_energy)...
                                        *ENVB.NomData.Bus.energy_trajectory(i);
    CurDist_B2{i}.Nomdata.Var_power = sum(DesignB.Var_energy_flow(1:ENVB.RealData.Battery.num,1:ENVB.RealData.Battery.num),2)...
                                        /(-DesignB.Maximum_output_energy)...
                                        *ENVB.NomData.Bus.energy_trajectory(i);
    CurDist_B2{i}.Nomdata.Sum_power = CurDist_B2{i}.Nomdata.Avg_power + CurDist_B2{i}.Nomdata.Var_power;
    
    Exp_B2{i}.Nomdata.Avg_Power = DesignB.Avg_energy_flow(1:ENVB.RealData.Battery.num,1:ENVB.RealData.Battery.num)...
                                        /(-DesignB.Maximum_output_energy)...
                                        *ENVB.NomData.Bus.energy_trajectory(i);
    Exp_B2{i}.Nomdata.Var_Power = DesignB.Var_energy_flow(1:ENVB.RealData.Battery.num,1:ENVB.RealData.Battery.num)...
                                        /(-DesignB.Maximum_output_energy)...
                                        *ENVB.NomData.Bus.energy_trajectory(i);
                                         
    for j = 1:ENVB.RealData.Battery.num 
        CurDist_B2{i}.Nomdata.Avg_Current(j) = CurDist_B2{i}.Nomdata.Avg_power(j)./CL.Bat{j}.volt;
        CurDist_B2{i}.Nomdata.Var_Current(j) = CurDist_B2{i}.Nomdata.Var_power(j)./CL.Bat{j}.volt;
        CurDist_B2{i}.Nomdata.Sum_Current(j) = CurDist_B2{i}.Nomdata.Sum_power(j)./CL.Bat{j}.volt;
    end
    CurDist_B2{i}.Nomdata.Bus_Current = CurDist_B2{i}.Nomdata.Bus_power./ENVB.NomData.Bus.volt;
end
for i = 1:length(ENVB.NomData.Bus.energy_trajectory)
    CurDist_B2{i}.Realdata.Avg_Current = CurDist_B2{i}.Nomdata.Avg_Current*ENVB.Base.curlim;
    CurDist_B2{i}.Realdata.Var_Current = CurDist_B2{i}.Nomdata.Var_Current*ENVB.Base.curlim;
    CurDist_B2{i}.Realdata.Sum_Current = CurDist_B2{i}.Nomdata.Sum_Current*ENVB.Base.curlim;
    CurDist_B2{i}.Realdata.Bus_Current = CurDist_B2{i}.Nomdata.Bus_Current*ENVB.Base.curlim;
    
    Exp_B2{i}.Realdata.Avg_Power = Exp_B2{i}.Nomdata.Avg_Power*ENVB.Base.power;
    Exp_B2{i}.Realdata.Var_Power = Exp_B2{i}.Nomdata.Var_Power*ENVB.Base.power;
    
    CurDist_B2{i}.Alpha.Avg = CurDist_B2{i}.Nomdata.Avg_Current/CurDist_B2{i}.Nomdata.Bus_Current;
    CurDist_B2{i}.Alpha.Var = CurDist_B2{i}.Nomdata.Var_Current/CurDist_B2{i}.Nomdata.Bus_Current;
    CurDist_B2{i}.Alpha.Sum = CurDist_B2{i}.Nomdata.Sum_Current/CurDist_B2{i}.Nomdata.Bus_Current;
end
save('CurDist_B2','CurDist_B2');
% Hardware experiments
save('Exp_B2','Exp_B2');
