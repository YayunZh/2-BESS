% %%  User Input
% global CL ENV;
% %User Input begin
% % please always input Baterry first, then capacitor and at last bus
% clear all;
% clc;
% % create 9 batteries, 1 bus and 1+1 cap
% for i = 1:9
%     Input{i} =  struct('sn','bat_1','ind',i,'volt',1,'ohmres',0.01,...
%        'curlim',1,'curlim_var',0.2,'curlim_mu',1,... % current lim,hete
%        'res',0.001,'res_mu',0.003,'res_var',0.001/3,...
%        'qlim',1,'qlim_mu',1,'qlim_var',0.2,...% electric charge
%        'eloss_coeff',0.12,'eloss_coeff_mu',0.12,'eloss_coeff_var',0.012);  % propotional for power process
% end
% % no power in cap
% Input{end+1} = struct('sn','bus1','ind',1,'volt',1,'curlim',-1); % index 10
% 
% Input{end+1} = struct('sn','cap1','ind',1,'volt',1,'capci',10000,'ohmres',0.001,'curlim',100);% index 11
% 
% Input{end+1} = struct('sn','cap2','ind',3,'volt',1,'capci',10000,'ohmres',0.001,'curlim',100);% index 13
% 
% % environment
% % simulation parameters
% ENV = struct('Constraint',[],'Sweep',[],'Fir_Conv',[]); %,'Sec_Conv',[]);   % Initialize the environment
% ENV.Constraint.lambda = 0.3;
% ENV.Sweep.Stat.Conv = 5;    % number of sec converter power limit to sweep
% ENV.Sweep.Stat.Bat = 4;     % number of battery types to sweep 
% ENV.Sweep.Conv.p_lim = linspace(0.001,0.2,ENV.Sweep.Stat.Conv);% 
% ENV.Sweep.Conv.e_lim = linspace(0.005,0.05,ENV.Sweep.Stat.Conv);%
% % layer 1
% ENV.Fir_Conv.Num = 9;       % the upper num lim of first-layer converter
% ENV.Fir_Conv.p_lim_single = 100;   % the upper power lim of first-layer converter, all var converters have the same power limit
% ENV.Fir_Conv.e_lim_single = 100;   % the upper power lim of first-layer converter, all var converters have the same power limit
% ENV.Fir_Conv.MC_trial = 2;  % 
% % % layer 2
% % ENV.Sec_Conv.Num = 9;       % the upper num lim of sec-layer converter
% % ENV.Sec_Conv.p_lim_single = 100;   % the upper power lim of sec-layer converter, all var converters have the same power limit
% % ENV.Sec_Conv.e_lim_single = 100;   % the upper power lim of sec-layer converter, all var converters have the same power limit
% % ENV.Sec_Conv.MC_trial = 2;  % takes 1 hours
% % Categorize the Component lists into Batteries, Cpacitors and Buses
% CL = struct('Bat',[],'Cap',[],'Bus',[],'Conn',[],'Stat',[]);   % Initialize the compoennt list struct
% for i = 1:size(Input,2)
%     Input{i}.ind = i;  % updated order is bat, bus and cap
%     if strcmp(Input{i}.sn(1:3),'bat')
%         CL.Bat{end+1} = Input{i};
%     elseif strcmp(Input{i}.sn(1:3),'bus')
%         CL.Bus{end+1} = Input{i};
%     elseif strcmp(Input{i}.sn(1:3),'cap')
%         CL.Cap{end+1} = Input{i};
%     else
%         fprintf('Error in input\n');
%         pause;
%     end
% end
% CL.Stat.Bat_num = size(CL.Bat,2);
% CL.Stat.Cap_num = size(CL.Cap,2);
% CL.Stat.Bus_num = size(CL.Bus,2);
% CL.Stat.Node_num = size(CL.Bat,2) + size(CL.Cap,2) + size(CL.Bus,2);
% for i = 1:CL.Stat.Cap_num   %Initialize the capactior votlages
%     CL.Cap{i}.volt = Inf;  
% end
% % Report the Component lists
% fprintf('%d%s\n%d%s\n%d%s\n',size(CL.Bat,2),' Batteries',size(CL.Cap,2),' SuperCapacitors',size(CL.Bus,2),' Buses');
% % used in sweep, need checking
% for i = 1:CL.Stat.Bat_num
%     ENV.Sweep.Bat{i}.curlim_mu = linspace(2,2,ENV.Sweep.Stat.Bat); % number of battery types to sweep 
%     ENV.Sweep.Bat{i}.curlim_var = linspace(0.1,0.4,ENV.Sweep.Stat.Bat);
%     ENV.Sweep.Bat{i}.res_mu = linspace(0,0,ENV.Sweep.Stat.Bat);
%     ENV.Sweep.Bat{i}.res_var = linspace(0,0,ENV.Sweep.Stat.Bat);
% end
% 
% P_direct_mat_in = zeros(CL.Stat.Bat_num + CL.Stat.Bus_num + CL.Stat.Cap_num - 1 , CL.Stat.Bat_num + CL.Stat.Bus_num + CL.Stat.Cap_num - 1 );  % 12 * 12 initialize direct power delivery interconnection matrix for P
% % Input the direct power delivery connection
% P_direct_mat_in(CL.Bat{1}.ind, CL.Bus{1}.ind) = 1;
% P_direct_mat_in(CL.Bat{2}.ind, CL.Bus{1}.ind) = 1;
% P_direct_mat_in(CL.Bat{3}.ind, CL.Bus{1}.ind) = 1;
% P_direct_mat_in(CL.Bat{4}.ind, CL.Bus{1}.ind) = 1;
% P_direct_mat_in(CL.Bat{5}.ind, CL.Bus{1}.ind) = 1;
% P_direct_mat_in(CL.Bat{6}.ind, CL.Bus{1}.ind) = 1;
% P_direct_mat_in(CL.Bat{7}.ind, CL.Bus{1}.ind) = 1;
% P_direct_mat_in(CL.Bat{8}.ind, CL.Bus{1}.ind) = 1;
% P_direct_mat_in(CL.Bat{9}.ind, CL.Bus{1}.ind) = 1;
% %P_direct_mat_in(CL.Cap{1}.ind, CL.Bus{1}.ind) = 1;
% % Input the diff power delivery connection
% % if i- j has conn, p(i,j) = 1
% P_diff_1_mat_in = zeros(CL.Stat.Bat_num + CL.Stat.Bus_num + CL.Stat.Cap_num - 1 , CL.Stat.Bat_num + CL.Stat.Bus_num + CL.Stat.Cap_num - 1 );
% P_diff_1_mat_in(CL.Bat{1}.ind, CL.Cap{1}.ind) = 1;
% P_diff_1_mat_in(CL.Bat{2}.ind, CL.Cap{1}.ind) = 1;
% P_diff_1_mat_in(CL.Bat{3}.ind, CL.Cap{1}.ind) = 1;
% P_diff_1_mat_in(CL.Bat{4}.ind, CL.Cap{1}.ind) = 1;
% P_diff_1_mat_in(CL.Bat{5}.ind, CL.Cap{1}.ind) = 1;
% P_diff_1_mat_in(CL.Bat{6}.ind, CL.Cap{1}.ind) = 1;
% P_diff_1_mat_in(CL.Bat{7}.ind, CL.Cap{1}.ind) = 1;
% P_diff_1_mat_in(CL.Bat{8}.ind, CL.Cap{1}.ind) = 1;
% P_diff_1_mat_in(CL.Bat{9}.ind, CL.Cap{1}.ind) = 1;
% 
% % P_diff_2_mat_in = zeros(CL.Stat.Bat_num + CL.Stat.Bus_num + CL.Stat.Cap_num - 1 , CL.Stat.Bat_num + CL.Stat.Bus_num + CL.Stat.Cap_num - 1 );
% % P_diff_2_mat_in(CL.Bat{1}.ind, CL.Cap{2}.ind) = 1;
% % P_diff_2_mat_in(CL.Bat{2}.ind, CL.Cap{2}.ind) = 1;
% % P_diff_2_mat_in(CL.Bat{3}.ind, CL.Cap{2}.ind) = 1;
% % P_diff_2_mat_in(CL.Bat{4}.ind, CL.Cap{2}.ind) = 1;
% % P_diff_2_mat_in(CL.Bat{5}.ind, CL.Cap{2}.ind) = 1;
% % P_diff_2_mat_in(CL.Bat{6}.ind, CL.Cap{2}.ind) = 1;
% % P_diff_2_mat_in(CL.Bat{7}.ind, CL.Cap{2}.ind) = 1;
% % P_diff_2_mat_in(CL.Bat{8}.ind, CL.Cap{2}.ind) = 1;
% % P_diff_2_mat_in(CL.Bat{9}.ind, CL.Cap{2}.ind) = 1;
% %User Input end
% 
% % Quick Check
% % InputCheck(CL,ENV);
% CL.Conn.direct = P_direct_mat_in;
% CL.Conn.diff_1 = P_diff_1_mat_in;
% % CL.Conn.diff_2 = P_diff_2_mat_in;
% % CL.Stat.VarConn_num = sum(sum(P_diff_2_mat_in));% 9
% 
% ENV.FakeBatt.branch_cur = ones(1,CL.Stat.Bat_num);% 1*9 [1]
% 
% save('CL_LSV_L1','CL');
% save('ENV_LSV_L1','ENV');

% % Task 5: Energy Design Converters 
% clear;
% load('ENV_LSV');                    
% load('CL_LSV');
% ENV5 = ENV;
% ENV5.Sweep.Stat.Bat = 8;       % number of battery types to sweep 
% ENV5.Sweep.Bat_energy_sum = 0;
% for i = 1:CL.Stat.Bat_num
%     ENV5.Sweep.Bat{i}.qlim_mu = linspace(1,1,ENV5.Sweep.Stat.Bat);
%     ENV5.Sweep.Bat{i}.qlim_var = [0.05, 0.1, 0.15, 0.2, 0.25, 0.5,0.75, 1]; %Capacity Variation
%     ENV5.Sweep.Bat_energy_sum = ENV5.Sweep.Bat_energy_sum + ENV5.Sweep.Bat{i}.qlim_mu(1);
% end
% ENV5.Sweep.Stat.Conv = 20;  % number of converter power limits to sweep

% load('ENV4_1');
% ENV5.Sweep.Conv.e_lim = ENV4_1.Sweep.Conv_energy_sum/ENV5.Fir_Conv.Num;
% 
% 
% % ENV5.Sweep.Conv.e_lim = linspace(0,1,ENV5.Sweep.Stat.Conv);
% ENV5.Sweep.Stat.Conv = length(ENV5.Sweep.Conv.e_lim);
% 
% ENV5.Fir_Conv.MC_trial = 20;  % (originally 100)
% tic
% ExpRes5 = func_lsv_converter_energy_design(CL,ENV5); %k0,k1,k2,:
% toc
% 
% save('ExpRes5','ExpRes5');
% save('ENV5','ENV5');

% %  % % Plot Task 5
% % clear;
% load('ENV4_trad');
% load('ExpRes4_trad');
% load('ENV5'); % Need to run task 4-1 ahead
% load('CL_LSV');
% load('ExpRes5');
% % func_plot_task_lsv(CL,ENV5,ExpRes5);
% % % %  Plot Task 5-1 Compare to LSH-PPP
% % load('ENV4_1');
% % load('ExpRes4');
% func_plot_task_lsv5_1(CL,ENV5,ExpRes5,ENV4_trad,ExpRes4_trad);
% % % %  Plot Task 5-2 Compare to LSH-PPP and FPP
% % load('ENV4_fpp','ENV4_fpp');
% % load('ExpRes4_fpp','ExpRes4_fpp');
% % func_plot_task_lsv_5_2(CL,ENV5,ExpRes5,ENV4_1,ExpRes4,ENV4_fpp,ExpRes4_fpp);
