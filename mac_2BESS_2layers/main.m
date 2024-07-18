clear;
%% Task 1 modulation of signal
% Parameters
ADC.f = 1;                % frequency
ADC.V_ref = 216;          % output voltage fundamental component peak
ADC.n_bat = 9;            % the number of batteries
ADC.n_v_level = 2 * ADC.n_bat + 1; % the number of voltage levels
ADC.n_t = 4 * ADC.n_bat + 2;  % the number of time intervals

ADC.num_samples = 1000;   % Number of samples to convert

battery_voltage = 24;

% Generate analog input signal
ADC.t = linspace(0, 1/ADC.f, ADC.num_samples + 1);
ADC.dt = 1/(ADC.f*ADC.num_samples);
analog_signal = ADC.V_ref * 1.01 * sin(2*pi*ADC.f*ADC.t);

% Simulate quantization (ADC conversion)
digital_signal = quantize(analog_signal,battery_voltage);
% Plot the analog and digital signals
plot(analog_signal, 'Color',[0.000000, 0.466667, 0.713725], 'LineWidth', 1.5);
hold on;
stairs(digital_signal, 'Color',[0.898039, 0.219608, 0.231373], 'LineWidth', 1);
hold off;
xlabel('Sample');
ylabel('Amplitude');
title('ADC Simulation');
legend('Analog Signal', 'Digital Signal');


% produce output and delta_t vector
delta_t = zeros(1,ADC.n_t); % Time corresponding to different voltage levels
output = zeros(1,ADC.n_t);
OUT.V_level = battery_voltage;

for i = 1: ceil(ADC.n_t/4)
    % OUT.value = round((i - 1) * battery_voltage,4);
    OUT.value_level = (i - 1) * OUT.V_level;

    indices = find(digital_signal(1:ceil(length(digital_signal)/4)) == OUT.value_level);

    if i ~= ceil(ADC.n_t/4)
        delta_t(i) = ADC.dt * (max(indices) - min(indices)+1);
        delta_t(2 * ADC.n_bat + 2 - i) = delta_t(i);
        delta_t(2 * ADC.n_bat + 1 + i) = delta_t(i);
        delta_t(4 * ADC.n_bat + 3 - i) = delta_t(i);
        output(i) = OUT.value_level;
        output(2 * ADC.n_bat + 2 - i) = output(i);
        output(2 * ADC.n_bat + 1 + i) = -output(i);
        output(4 * ADC.n_bat + 3 - i) = -output(i);
    else
        delta_t(i) = 2 * ADC.dt * (max(indices) - min(indices));
        delta_t(2 * ADC.n_bat + 1 + i) = delta_t(i);
        output(i) = OUT.value_level;
        output(2 * ADC.n_bat + 1 + i) = -output(i);
    end
end

delta_t_level = delta_t / sum(delta_t); % proportion

save('time_intervals','delta_t');
save('voltage_levels','output');
save('ADC','ADC');

%% Task 2 User Input
clear;
load('voltage_levels');
load('time_intervals');
load('ADC');

global SYS TEST;

% create batteries, HB, 1 bus 
for i = 1:ADC.n_bat
    Input{i} =  struct('sn','bat_1','ind',i,'volt',24,'ohmres',0.01,...
       'curlim',25,'curlim_var',2.5,'curlim_mu',25,... % current lim,heterogeneity
       'capalim',250,'capalim_mu',250,'capalim_var',0.1);% battery capacity
end
Input{end+1} = struct('sn','bus1','ind',1,'volt',1,'curlim',-1); 

%%
% % TEST environment
% % simulation parameters
TEST = struct('Constraint',[],'Sweep',[],'Conv',[]);   % Initialize the TESTironment
% TEST.Constraint.lambda = 0.3;
TEST.Sweep.Stat.Conv = 10;   % number of var converter power limit to sweep
TEST.Sweep.Stat.Bat = 5;     % number of battery types to sweep 
TEST.Sweep.Conv.p_lim = linspace(0.1,1,TEST.Sweep.Stat.Conv);%linspace(0.01,0.9,TEST.Sweep.Stat.Conv);% 

TEST.Conv.Num = ADC.n_bat - 1;   % the upper num lim of dense converter
TEST.Conv.p_lim_singlevar = 1;   % the upper power lim of converter, all dense converters have the same power limit
TEST.Conv.MC_trial = 100;  

% Categorize the Component lists into Batteries, Cpacitors and Buses
SYS = struct('Bat',[],'Cap',[],'Bus',[],'Conn',[],'Stat',[]);   % Initialize the compoennt list struct
for i = 1:size(Input,2)
    Input{i}.ind = i;  % updated order is bat, bus and cap
    if strcmp(Input{i}.sn(1:3),'bat')
        SYS.Bat{end+1} = Input{i};
    elseif strcmp(Input{i}.sn(1:3),'bus')
        SYS.Bus{end+1} = Input{i};
    elseif strcmp(Input{i}.sn(1:3),'cap')
        SYS.Cap{end+1} = Input{i};
    else
        fprintf('Error in input\n');
        pause;
    end
end

SYS.Stat.Bat_num = size(SYS.Bat,2);
SYS.Stat.Cap_num = size(SYS.Cap,2);
SYS.Stat.Bus_num = size(SYS.Bus,2);
SYS.Stat.Node_num = size(SYS.Bat,2);
SYS.Stat.Delta_t_num = size(delta_t,2);
SYS.Stat.Delta_t = delta_t;
SYS.Stat.Output = output;

% Report the Component lists
fprintf('%d%s\n%d%s\n%d%s\n',size(SYS.Bat,2),' Batteries',size(SYS.Cap,2),' SuperCapacitors',size(SYS.Bus,2),' Buses');

for i = 1:SYS.Stat.Bat_num
    TEST.Sweep.Bat{i}.curlim_mu = linspace(1*25,1*25,TEST.Sweep.Stat.Bat); % number of battery types to sweep 
    TEST.Sweep.Bat{i}.curlim_var = linspace(0.05*25,0.25*25,TEST.Sweep.Stat.Bat);
end

% if i- j has conn, p(i,j) = 1
P_diff_mat_in = zeros(SYS.Stat.Bat_num, SYS.Stat.Bat_num); % initialize differential power delivery interconnection matrix for P  
for i = 1:(SYS.Stat.Bat_num-1)
    P_diff_mat_in(SYS.Bat{i}.ind, SYS.Bat{i+1}.ind) = 1;
end
% User Input end
SYS.Conn.diff = P_diff_mat_in;

save('SYS','SYS');
save('TEST','TEST');


%% Task 3 Connection Design (layer 1 Converter & battery use)
clear;
load('TEST');
load('SYS');
TEST3 = TEST;
SYS3 = SYS;

TEST3.Conv.flow = 1; % chosen to have constant(1) or variable(0) power flow during an ac cycle

TEST3.Fir_Conv.Num = 3; % number of fist layer converter
TEST3.Conv.partition = 1; % how to partition converters to achieve the economic of scales
TEST3.Fir_Conv.Position = 0; % select experience positions to save time(1) or re-optimisation(0)

TEST3.Sweep.Stat.Conv = 5;   % number of var converter power limit to sweep
TEST3.Sweep.Conv.p_lim = linspace(0.1,1,TEST.Sweep.Stat.Conv);%linspace(0.01,0.9,TEST.Sweep.Stat.Conv);% 
TEST3.Sweep.Stat.Bat = 5;     % number of battery types to sweep - heterogeneity
for i = 1:SYS3.Stat.Bat_num
    TEST3.Sweep.Bat{i}.curlim_mu = linspace(1*25,1*25,TEST3.Sweep.Stat.Bat); 
    TEST3.Sweep.Bat{i}.curlim_var = linspace(0.05*25,0.25*25,TEST3.Sweep.Stat.Bat);
end

save('SYS3','SYS3');
save('TEST3','TEST3');

tic
Opt_connection = func_connection_design(SYS3,TEST3);
toc

save('Opt_connection.mat','Opt_connection');

%% Task 4 MC simulation
clear;
load('TEST3');
load('SYS3');
TEST4 = TEST3;
SYS4 = SYS3;

TEST4.Conv.MC_trial = 2; % number of MC simulation

save('SYS4','SYS4');
save('TEST4','TEST4');
tic
Opt_MC = func_MC(SYS4,TEST4);
toc

save('MC.mat','Opt_MC');


% %% Task 5 Output PQ
% clear;
% load('TEST3');
% load('SYS3');
% TEST5 = TEST3;
% SYS5 = SYS3;
% 
% TEST5.Sweep.Stat.Conv = 10;   % number of var converter power limit to sweep
% TEST5.Sweep.Conv.p_lim = linspace(0.1,1,TEST5.Sweep.Stat.Conv);%linspace(0.01,0.9,TEST.Sweep.Stat.Conv);% 
% TEST5.Sweep.Stat.Bat = 5;     % number of battery types to sweep - heterogeneity
% for i = 1:SYS3.Stat.Bat_num
%     TEST5.Sweep.Bat{i}.curlim_mu = linspace(1*25,1*25,TEST5.Sweep.Stat.Bat); 
%     TEST5.Sweep.Bat{i}.curlim_var = linspace(0.05*25,0.25*25,TEST5.Sweep.Stat.Bat);
% end
% 
% TEST5.Sweep.Stat.phi = 100; % number of phase shift
% TEST5.Sweep.phi = linspace(-0.99,0,TEST5.Sweep.Stat.phi);
% TEST5.t_delay = 1; 
% 
% TEST5.Fir_Conv.Position = 1;
% 
% save('SYS5','SYS5');
% save('TEST5','TEST5');
% tic
% Opt_PQ = func_PQ(SYS5,TEST5);
% toc

%%
function digital_signal = quantize(analog_signal, battery_voltage)
    % Calculate LSB (Least Significant Bit)
    v_level = battery_voltage;
    temp_signal = fix(analog_signal / v_level) * v_level;
    correction = 1;% / max(temp_signal);
    digital_signal = temp_signal * correction;   
end