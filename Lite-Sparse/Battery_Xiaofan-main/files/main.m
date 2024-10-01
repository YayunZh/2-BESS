% clear;
% load('CL');
% CL = Flatten_battery_stat(CL);

%%  User Input
global CL ENV;
%User Input begin
% please always input Baterry first, then capacitor and at last bus
clear all;
clc;
% create 9 batteries, 1 bus and 1 cap
for i = 1:9
    Input{i} =  struct('sn','bat_1','ind',i,'volt',1,'ohmres',0.01,...
       'curlim',1,'curlim_var',0.2,'curlim_mu',1,... % current lim,hete
       'res',0.001,'res_mu',0.003,'res_var',0.001/3,...
       'qlim',1,'qlim_mu',1,'qlim_var',0.2,...% electric charge
       'eloss_coeff',0.12,'eloss_coeff_mu',0.12,'eloss_coeff_var',0.012);  % propotional for power process
end
% no power in cap
Input{end+1} = struct('sn','bus1','ind',1,'volt',1,'curlim',-1); % index 10
Input{end+1} = struct('sn','cap1','ind',1,'volt',1,'capci',10000,'ohmres',0.001,'curlim',100);% index 11

% environment
% simulation parameters
ENV = struct('Constraint',[],'Sweep',[],'Avg_Conv',[],'Var_Conv',[]);   % Initialize the environment
ENV.Constraint.lambda = 0.3;
ENV.Sweep.Stat.Conv = 5;    % number of var converter power limit to sweep
ENV.Sweep.Stat.Bat = 4;     % number of battery types to sweep 
ENV.Sweep.Conv.p_lim = linspace(0.001,0.2,ENV.Sweep.Stat.Conv);% 
ENV.Sweep.Conv.e_lim = linspace(0.005,0.05,ENV.Sweep.Stat.Conv);%

% layer 1
ENV.Avg_Conv.Num = 3;      % the upper num lim of averaging converter
ENV.Avg_Conv.ActNum = 3;   % the actural num of averaging converter can be <= the num lim of averaging converter, becuase some of them does not deliver power
ENV.Avg_Conv.p_lim_vec = 100*ones(1,ENV.Avg_Conv.Num);  % the upper power lim of averaging converter
ENV.Avg_Conv.e_lim_vec = 100*ones(1,ENV.Avg_Conv.Num);  % the upper energy lim of averaging converter

% layer 2
ENV.Var_Conv.Num = 8;       % the upper num lim of var converter
ENV.Var_Conv.p_lim_singlevar = 100;   % the upper power lim of var converter, all var converters have the same power limit
ENV.Var_Conv.e_lim_singlevar = 100;   % the upper power lim of var converter, all var converters have the same power limit
ENV.Var_Conv.MC_trial = 2;  % takes 1 hours

%Categorize the Component lists into Batteries, Capacitors and Buses
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
CL.Stat.Node_num = size(CL.Bat,2) + size(CL.Cap,2) + size(CL.Bus,2);
for i = 1:CL.Stat.Cap_num   %Initialize the capactior votlages
    CL.Cap{i}.volt = Inf;  
end
% Report the Component lists
fprintf('%d%s\n%d%s\n%d%s\n',size(CL.Bat,2),' Batteries',size(CL.Cap,2),' SuperCapacitors',size(CL.Bus,2),' Buses');
for i = 1:CL.Stat.Bat_num
    ENV.Sweep.Bat{i}.curlim_mu = linspace(2,2,ENV.Sweep.Stat.Bat); % number of battery types to sweep 
    ENV.Sweep.Bat{i}.curlim_var = linspace(0.1,0.4,ENV.Sweep.Stat.Bat);
    ENV.Sweep.Bat{i}.res_mu = linspace(0,0,ENV.Sweep.Stat.Bat);
    ENV.Sweep.Bat{i}.res_var = linspace(0,0,ENV.Sweep.Stat.Bat);
end


P_direct_mat_in = zeros(CL.Stat.Bat_num + CL.Stat.Bus_num , CL.Stat.Bat_num + CL.Stat.Bus_num );  % initialize direct power delivery interconnection matrix for P
P_diff_mat_in = zeros(CL.Stat.Bat_num + CL.Stat.Bus_num, CL.Stat.Bat_num + CL.Stat.Bus_num);    % initialize differential power delivery interconnection matrix for P  
% Input the direct power delivery connection
% B1-B9,C1 with Bus 10*10
P_direct_mat_in(CL.Bat{1}.ind, CL.Bus{1}.ind) = 1;
P_direct_mat_in(CL.Bat{2}.ind, CL.Bus{1}.ind) = 1;
P_direct_mat_in(CL.Bat{3}.ind, CL.Bus{1}.ind) = 1;
P_direct_mat_in(CL.Bat{4}.ind, CL.Bus{1}.ind) = 1;
P_direct_mat_in(CL.Bat{5}.ind, CL.Bus{1}.ind) = 1;
P_direct_mat_in(CL.Bat{6}.ind, CL.Bus{1}.ind) = 1;
P_direct_mat_in(CL.Bat{7}.ind, CL.Bus{1}.ind) = 1;
P_direct_mat_in(CL.Bat{8}.ind, CL.Bus{1}.ind) = 1;
P_direct_mat_in(CL.Bat{9}.ind, CL.Bus{1}.ind) = 1;
%P_direct_mat_in(CL.Cap{1}.ind, CL.Bus{1}.ind) = 1;
% Input the diff power delivery connection
% if i- j has conn, p(i,j) = 1
P_diff_var_mat_in = zeros(CL.Stat.Bat_num + CL.Stat.Bus_num, CL.Stat.Bat_num + CL.Stat.Bus_num);
P_diff_var_mat_in(CL.Bat{1}.ind, CL.Bat{2}.ind) = 1;
P_diff_var_mat_in(CL.Bat{2}.ind, CL.Bat{3}.ind) = 1;
P_diff_var_mat_in(CL.Bat{3}.ind, CL.Bat{4}.ind) = 1;
P_diff_var_mat_in(CL.Bat{4}.ind, CL.Bat{5}.ind) = 1;
P_diff_var_mat_in(CL.Bat{5}.ind, CL.Bat{6}.ind) = 1;
P_diff_var_mat_in(CL.Bat{6}.ind, CL.Bat{7}.ind) = 1;
P_diff_var_mat_in(CL.Bat{7}.ind, CL.Bat{8}.ind) = 1;
P_diff_var_mat_in(CL.Bat{8}.ind, CL.Bat{9}.ind) = 1;
%P_diff_var_mat_in(CL.Bat{9}.ind, CL.Cap{1}.ind) = 1;
% the traditional DPP architectures
% User Input end

% Quick Check
InputCheck(CL,ENV);
CL.Conn.direct = P_direct_mat_in;
CL.Conn.diff_var = P_diff_var_mat_in;
CL.Conn.diff_avg = zeros(CL.Stat.Bat_num + CL.Stat.Bus_num, CL.Stat.Bat_num + CL.Stat.Bus_num);
CL.Stat.VarConn_num = sum(sum(P_diff_var_mat_in));% 9

ENV.FakeBatt.branch_cur = ones(1,CL.Stat.Bat_num);% 1*9 [1]

save('CL','CL');
save('ENV','ENV');

%% Task 4 Energy Design Layer 1 Converters
clear;
load('ENV');
load('CL');
ENV4 = ENV;
ENV4.Avg_Conv.trial_num = 2;%10000;  % the numer of random searches to find the optimal placement of averaging converter
ENV4.Avg_Conv.partition = 1; % how to partition the layer 1 converters to achieve the economic of scales
% partition the averaging converter into groups  
% Energy Design Average Converters
% if (isfile('OptAvg_7e.mat'))&&(isfile('OptAvg_5e.mat'))&&(isfile('OptAvg_3e.mat'))
if (isfile('OptAvg_3e.mat'))
    u_in = input('To Save Time, We can Use the Previous 3/5/7 Avg Conv for Var Converter Design, type 1(Y)/0(N): ');
    if (u_in == 1)
%         load('OptAvg_7.mat'); 
%         load('OptAvg_5.mat'); 
        load('OptAvg_3e.mat');
    else
        fprintf('No Previous Ave Conv Design, Need Spend Time to Design Avg Conv, Check the ENV to make sure the setting\n');
        OptAvg_energy = func_average_converter_energy_design(CL,ENV4);
    end
else
    tic
    OptAvg_energy = func_average_converter_energy_design(CL,ENV4);
    toc
end
% postprocessing: sort the converter rating in ascending orders
[OptAvg_energy.Conv_energy_rating_partition, OptAvg_energy.Conv_energy_rating_partition_mat] ...
    = func_rating_partition(abs(OptAvg_energy.Avg_energy_flow), ENV4.Avg_Conv.partition); 
% postprocessing: update the primary converter design results to ENV4
ENV4.Avg_Conv.e_lim_vec = OptAvg_energy.Conv_energy_rating_partition;
ENV4.Avg_Conv.e_lim_mat = OptAvg_energy.Conv_energy_rating_partition_mat;
ENV4.Avg_Conv.ActNum = sum(sum(OptAvg_energy.Conv_energy_rating_partition_mat > 0))/2;

fprintf('The battery utilization is\n');    % ~75%
ENV4.Avg_Conv.bat_uratio = abs(OptAvg_energy.output_energy)/sum(OptAvg_energy.Bat_energy_rating);
ENV4.Avg_Conv.bat_uratio
fprintf('The Energy Utilization ratio of the Average Converters\n');  % ~100%
ENV4.Avg_Conv.energy_uratio = sum(sum(abs(OptAvg_energy.Avg_energy_flow)))/sum(sum(ENV4.Avg_Conv.e_lim_mat));
ENV4.Avg_Conv.energy_uratio
% for i = 1:CL.Stat.Bat_num
%     ENV4.Sweep.Bat{i}.eloss_coeff_mu = linspace(0.12,0.12,ENV4.Sweep.Stat.Bat);
%     ENV4.Sweep.Bat{i}.eloss_coeff_var = linspace(0.002,0.012,ENV4.Sweep.Stat.Bat);
% end
save('OptAvg_3e','OptAvg_energy');
save('ENV4','ENV4');

% % % Task 4-1: Energy Design Layer 2 Converters 
% % Need to run task4 ahead
clear;
load('ENV4');                    
load('OptAvg_3e');
load('CL');
ENV4_1 = ENV4;
ENV4_1.Sweep.Stat.Bat = 8;       % number of battery types to sweep 
ENV4_1.Sweep.Bat_energy_sum = 0;
for i = 1:CL.Stat.Bat_num
    ENV4_1.Sweep.Bat{i}.qlim_mu = linspace(1,1,ENV4_1.Sweep.Stat.Bat);
    ENV4_1.Sweep.Bat{i}.qlim_var = [0.05, 0.1, 0.15, 0.2, 0.25, 0.5,0.75, 1]; %Capacity Variation
    ENV4_1.Sweep.Bat_energy_sum = ENV4_1.Sweep.Bat_energy_sum + ENV4_1.Sweep.Bat{i}.qlim_mu(1);
end
ENV4_1.Sweep.Stat.Conv = 20;  % number of var converter power limits to sweep
%ENV4_1.Sweep.Stat.Conv = 3;  % number of var converter power limits to sweep
ENV4_1.Sweep.Conv_energy_sum_primary = sum(ENV4_1.Avg_Conv.e_lim_vec);% 3 avg conv
secondary_energy_sum_1 = linspace(0,2,10)*sum(ENV4_1.Avg_Conv.e_lim_vec);
secondary_energy_sum_2 =  - sum(ENV4_1.Avg_Conv.e_lim_vec);
temp = sort(unique([secondary_energy_sum_1,secondary_energy_sum_2]));
ENV4_1.Sweep.Conv_energy_sum_secondary = [];
for i = 1:length(temp)
   if (temp(i) >= 0)
       ENV4_1.Sweep.Conv_energy_sum_secondary(end+1) = temp(i);
   end
end
ENV4_1.Sweep.Conv_energy_sum = ENV4_1.Sweep.Conv_energy_sum_primary + ENV4_1.Sweep.Conv_energy_sum_secondary; % 1*20
ENV4_1.Sweep.Stat.Conv = length(ENV4_1.Sweep.Conv_energy_sum);
% var 
ENV4_1.Sweep.Conv.e_lim = (ENV4_1.Sweep.Conv_energy_sum - ENV4_1.Sweep.Conv_energy_sum_primary)/ENV4_1.Var_Conv.Num;
ENV4_1.Var_Conv.MC_trial = 20;  % (originally 100)
% SelfCheck_ENV4_1(CL,ENV4_1);
tic
sort = 1;
ExpRes4 = func_var_converter_energy_design(CL,ENV4_1,sort); %k0,k1,k2,: 1-sort, 0:unsort
toc

if(sort)
    save('ExpRes4_sort','ExpRes4');
    save('ENV4_1_sort','ENV4_1');
else
    save('ExpRes4_unsort','ExpRes4');
    save('ENV4_1_unsort','ENV4_1');
end

%  Plot Task 4-1
clear;
load('ENV4_1_sort'); % Need to run task 4-1 ahead
load('CL');
load('ExpRes4_sort');
func_plot_task4_1(CL,ENV4_1,ExpRes4);

% %% Task 4-4: Compare to C-PPP
% clear;
% load('ENV4_1'); % Need to run task 4-1 ahead
% load('CL');
% load('ExpRes4');
% ENV4_trad = ENV4_1;
% ENV4_trad.Avg_Conv.e_lim_vec = zeros(size(ENV4_trad.Avg_Conv.e_lim_vec));
% ENV4_trad.Avg_Conv.e_lim_mat = zeros(size(ENV4_trad.Avg_Conv.e_lim_mat));
% ENV4_trad.Sweep.Conv.e_lim = ENV4_trad.Sweep.Conv_energy_sum/ENV4_trad.Var_Conv.Num;
% ENV4_trad.Var_Conv.MC_trial = 20;  % takes 1 hours (originally 100)
% % run('MC_2layer.m');
% % OptAvg_energy_trad = OptAvg_energy;
% % OptAvg_energy_trad.Avg_energy_flow = zeros(CL.Stat.Node_num,CL.Stat.Node_num);
% % OptAvg_energy_trad.outputenergy = 0;
% % OptAvg_energy_trad.Conv_energy_rate = 0;
% tic
% sort=1;
% ExpRes4_trad = func_var_converter_energy_design(CL,ENV4_trad,sort);
% toc
% save('ExpRes4_trad','ExpRes4_trad');
% save('ENV4_trad','ENV4_trad');

% %%  Plot Task 4-4
% clear;
% load('ENV4_1'); % Need to run task 4-1 ahead
% load('CL');
% load('ExpRes4');
% load('ExpRes4_trad');
% load('ENV4_trad');
% func_plot_task4_4(CL,ENV4_1,ExpRes4,ENV4_trad,ExpRes4_trad);

% % % 
% % %% Task 4-5: Compare to C-PPP and FPP
% % clear;
% % load('CL');
% % load('ENV4_trad'); % Need to run task 4-1 ahead
% % load('ExpRes4_trad');
% % load('ENV4_1');
% % load('ExpRes4');
% ENV4_fpp = ENV4_trad;
% % ENV1_fpp.Avg_Conv.p_lim_vec = zeros(size(ENV1_fpp.Avg_Conv.p_lim_vec));
% % ENV1_fpp.Avg_Conv.p_lim_mat = zeros(size(ENV1_fpp.Avg_Conv.p_lim_mat));
% ENV4_fpp.Var_Conv.e_lim_singlevar = 0;
% ENV4_fpp.Sweep.Conv.e_lim = ENV4_fpp.Sweep.Conv_energy_sum/ENV4_fpp.Var_Conv.Num;
% ENV4_fpp.Var_Conv.MC_trial = 20;  % takes 1 hours (originally 12)
% tic
% ExpRes4_fpp = func_fpp_converter_energy_design(CL,ENV4_fpp);
% toc
% save('ENV4_fpp','ENV4_fpp');
% save('ExpRes4_fpp','ExpRes4_fpp');
% 
% %%  Plot Task 4-5
% clear;
% load('CL');
% load('ENV4_trad');    % Need to run task 4-4 ahead
% load('ExpRes4_trad');
% load('ENV4_1');
% load('ExpRes4');
% load('ENV4_fpp','ENV4_fpp');
% load('ExpRes4_fpp','ExpRes4_fpp');
% func_plot_task4_5(CL,ENV4_1,ExpRes4,ENV4_trad,ExpRes4_trad,ENV4_fpp,ExpRes4_fpp);
% %func_plot_task4_5_fakepower(CL,ENV4_1,ExpRes4,ENV4_trad,ExpRes4_trad,ENV4_fpp,ExpRes4_fpp);
