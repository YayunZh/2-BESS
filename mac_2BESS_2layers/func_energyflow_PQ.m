function OptRes = func_energyflow_PQ(SYS, TEST)
    
    P_Sec_mat_in = SYS.Conn.diff;
    P_Sec_mat = (P_Sec_mat_in > 0) | (P_Sec_mat_in > 0)';
    P_lim = TEST.Conv.p_lim_singlevar;
    Conv_Fir_num = TEST.Fir_Conv.Num; % First layer converter number

    Node_num = SYS.Stat.Bat_num;
    T_num = SYS.Stat.Delta_t_num;
    Delta_t = SYS.Stat.Delta_t;
    Output_waveform = SYS.Stat.Output;

    % the sum of the intrinsic power(current) capability of the batteries
    lim_sum = 0;
    for i = 1:Node_num
        lim_sum = lim_sum + SYS.Bat{i}.curlim;
    end

    Conv_lim = P_lim * lim_sum / (Node_num - 1 + Conv_Fir_num);

    yalmip('clear');
    % First optimization - battery power utilization
    ops = sdpsettings('solver', 'gurobi', ...    
                      'gurobi.Threads', 8, ...
                      'gurobi.MIPGap', 0.01, ... 
                      'gurobi.TimeLimit', 720, ...
                      'gurobi.Presolve', 2, ...
                      'gurobi.Cuts', 2,...
                      'gurobi.Heuristics', 1);

    I_max = sdpvar(1, 1, 'full');                 
    I_B = sdpvar(Node_num, T_num, 'full'); % Battery current matrix
    Q_B = sdpvar(Node_num, T_num, 'full'); % Battery charge matrix
    Q_C = sdpvar(Node_num, T_num, 'full'); % Converter charge matrix
    Q_L = sdpvar(Node_num, T_num, 'full'); % Converter charge matrix
    I_C_f = sdpvar(Node_num, Node_num, T_num, 'full'); % average Converter current matrix
    I_C_s = sdpvar(Node_num, Node_num, T_num, 'full'); % variance Converter current matrix
    Bat_use = binvar(Node_num, T_num, 'full'); % 
  
    if(TEST.Fir_Conv.Position == 0)
    P_Fir_mat = binvar(Node_num, Node_num, 'full');
    else
        P_Fir_mat = zeros(Node_num,Node_num);
        for i = 0:Conv_Fir_num-1
            P_Fir_mat(i+1,Node_num-i) = 1;
            P_Fir_mat(Node_num-i,i+1) = 1;
        end
    end

    Objective = sum(sum(Q_L));
    Constraints = [];


    %% constraints

    % [1] the constraint for charge matrix
    for i = 1:Node_num
        for j = 1:T_num
        Constraints = [Constraints, Q_B(i, j) == I_B(i, j) * Delta_t(j)];
        Constraints = [Constraints, Q_C(i, j) == (sum(I_C_s(i, :, j)) + sum(I_C_f(i, :, j))) * Delta_t(j)];
        Constraints = [Constraints, Q_L(i,j) == Q_C(i,j) + Q_B(i,j)];
        end
    end

    % [2] the constriant for output and phase shift
    Constraints = [Constraints, I_max >= 0];
    f = 1;
    M = 100;
    t_delay = TEST.t_delay; % phase shift

    t = linspace(0, 1/f, 1000);
    I_flow = I_max * sin(2 * pi * f * t); % assume sinusoidal wave
    t_extended = [t - 1/f, t, t + 1/f];
    I_flow_extended = [I_flow, I_flow, I_flow];

    for k = 1:T_num
        t_current(k) = sum(Delta_t(1:k)); 
        t_offset_start(k) = find(t_extended >= t_current(k) - Delta_t(k) - t_delay, 1);
        t_offset_end(k) = find(t_extended >= t_current(k) - t_delay, 1);
        T(k) = (t_offset_start(k)+t_offset_end(k))/2;
        I_values(k) = abs(I_flow_extended(fix((t_offset_start(k)+t_offset_end(k))/2))); % time k - corresponding current value

        for i = 1:Node_num
         Constraints = [Constraints, Q_L(i, k) - I_values(k) * Delta_t(k) <= M * (1 - Bat_use(i, k))];
         Constraints = [Constraints, I_values(k) * Delta_t(k) - Q_L(i, k) <= M * (1 - Bat_use(i, k))];
         end
    end

    % [3] the constriant for differential power processing links and Converter Power(current) Limit
   
    % first layer
    for i = 1:Node_num
        for j = i:Node_num
            for k = 1:T_num
            Constraints = [Constraints, I_C_f(i, j, k) <=   P_Fir_mat(i,j) * Conv_lim];
            Constraints = [Constraints, I_C_f(i, j, k) >= - P_Fir_mat(i,j) * Conv_lim];
            Constraints = [Constraints, I_C_f(i, j, k) == - I_C_f(j, i, k)];
            end
        end
    end

        % the constriant for first layer converter number
        if(TEST.Fir_Conv.Position == 0)
        for i = 1:Node_num
            for j = i:Node_num
            Constraints = [Constraints, P_Fir_mat(i,j) == P_Fir_mat(j,i)];
            end
        end
        Constraints = [Constraints, Conv_Fir_num == sum(sum(P_Fir_mat))/2];
        end

    % second layer
    for i = 1:Node_num
        for j = i:Node_num
            for k = 1:T_num
            Constraints = [Constraints, I_C_s(i, j, k) <=   P_Sec_mat(i,j) * Conv_lim];
            Constraints = [Constraints, I_C_s(i, j, k) >= - P_Sec_mat(i,j) * Conv_lim];
            Constraints = [Constraints, I_C_s(i, j, k) == - I_C_s(j, i, k)];
            end
        end
    end

    % [4] the constriant for battery Power(current) Limit
    for i = 1:Node_num
        for k = 1:T_num
        Constraints = [Constraints, SYS.Bat{i}.curlim >= I_B(i, k) - sum(I_C_f(i, :, k)) - sum(I_C_s(i, :, k))];
        Constraints = [Constraints, I_B(i, k) - sum(I_C_f(i, :, k)) - sum(I_C_s(i, :, k)) >= -SYS.Bat{i}.curlim];
        end
    end

    % [5] the constriant for SOC balancing
    for i = 2:Node_num
        Constraints = [Constraints, (sum(Q_B(i, :)) - sum(Q_C(i, :))) / (SYS.Bat{i}.capalim) == ...
                       (sum(Q_B(1, :)) - sum(Q_C(1, :))) / (SYS.Bat{1}.capalim)];
        SOC(1) = 1-((sum(Q_B(1, :)) - sum(Q_C(1, :))) / (SYS.Bat{1}.capalim));
        SOC(i) = 1-((sum(Q_B(i, :)) - sum(Q_C(i, :))) / (SYS.Bat{i}.capalim)) ;
    end

    % [6] the constriant for constant power flow of converters
    if (TEST.Conv.flow == 1)
        for i = 1:Node_num
            for j = i:Node_num
                for k = 1:T_num
                Constraints = [Constraints, I_C_s(i, j, k) == I_C_s(i, j, 1)];
                Constraints = [Constraints, I_C_f(i, j, k) == I_C_f(i, j, 1)];
                end
            end
        end
    end

    % [7] the constriant for output voltage and current
    for k = 1:T_num
        Constraints = [Constraints, sum(Bat_use(:, k)) == abs(floor(Output_waveform(k) / 24))];
    end

        % battery connect to the load
        M = 100;
        for k = 1:T_num
            for i = 1:Node_num
            Constraints = [Constraints, I_B(i, k) <= M * Bat_use(i, k)];
            Constraints = [Constraints, I_B(i, k) >= -M * Bat_use(i, k)];
            end
        end

    optimize(Constraints, -Objective, ops);

    OptRes.current_capability = lim_sum;
    OptRes.I_B = value(I_B);
    OptRes.I_max = value(I_max);
    OptRes.I_values = value(I_values);
    OptRes.I_C_s = value(I_C_s);
    OptRes.I_C_f = value(I_C_f);

    OptRes.Q_B = value(Q_B);
    OptRes.Q_C = value(Q_C);
    OptRes.Q_L = value(Q_L);
    OptRes.SOC = value(SOC);
    OptRes.Bat_use = value(Bat_use);
    OptRes.Fir_mat = value(P_Fir_mat);

    OptRes.T_start = value(t_offset_start);
    OptRes.T_end = value(t_offset_end);
    OptRes.T = value(T);
    OptRes.T_current = value(t_current);

end
