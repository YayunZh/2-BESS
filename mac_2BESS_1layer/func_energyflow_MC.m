function OptRes = func_energyflow_MC(Bat_Info_MC, SYS, TEST, QL_connection)
    
    P_mat_in = SYS.Conn.diff;
    P_mat = (P_mat_in > 0) | (P_mat_in > 0)';
    P_lim = TEST.Conv.p_lim_singlevar;
    
    Node_num = SYS.Stat.Bat_num;
    T_num = SYS.Stat.Delta_t_num;
    Delta_t = SYS.Stat.Delta_t;
    Output_waveform = SYS.Stat.Output;

    for i = 1: Node_num
        SYS.Bat{i}.curlim = Bat_Info_MC(i); % update batteries current limit 
        SYS.Bat{i}.capalim = Bat_Info_MC(i) * 10;
    end

    % the sum of the intrinsic power(current) capability of the batteries
    lim_sum = 0;
    for i = 1:Node_num
        lim_sum = lim_sum + SYS.Bat{i}.curlim;
    end

    Conv_lim = P_lim * lim_sum / (Node_num - 1);

    yalmip('clear');
    % First optimization - battery power utilization
    ops = sdpsettings('solver', 'gurobi', ...    
                      'gurobi.Threads', 8, ...
                      'gurobi.MIPGap', 0.01, ... 
                      'gurobi.TimeLimit', 720, ...
                      'gurobi.Presolve', 2, ...
                      'gurobi.Cuts', 2,...
                      'gurobi.Heuristics', 1);
                     
    I_B = sdpvar(Node_num, T_num, 'full'); % Battery current matrix
    I_L = sdpvar(T_num, 1, 'full');
    Q_B = sdpvar(Node_num, T_num, 'full'); % Battery charge matrix
    Q_C = sdpvar(Node_num, T_num, 'full'); % Converter charge matrix
    I_C = sdpvar(Node_num, Node_num, T_num, 'full'); % variance Converter current matrix
    Bat_use = QL_connection;
  
    Objective = sum(I_L);
    Constraints = [];


    %% constraints

    % [1] the constraint for charge matrix
    for i = 1:Node_num
        for j = 1:T_num
        Constraints = [Constraints, Q_B(i, j) == I_B(i, j) * Delta_t(j)];
        Constraints = [Constraints, Q_C(i, j) == sum(I_C(i, :, j)) * Delta_t(j)];
        end
    end

    % [2] the constriant for KCL 
    M = 100;
    for k = 1:T_num
         for i = 1:Node_num
         Constraints = [Constraints, sum(I_C(i, :, k)) + I_B(i,k) - I_L(k) <= M * (1 - Bat_use(i, k))];
         Constraints = [Constraints, I_L(k) - sum(I_C(i, :, k)) - I_B(i,k) <= M * (1 - Bat_use(i, k))];
         end
    end

    % [3] the constriant for differential power processing links and Converter Power(current) Limit

    if(TEST.Conv.partition==1)
        for i = 1:Node_num
            for j = i:Node_num
                for k = 1:T_num
                Constraints = [Constraints, I_C(i, j, k) <=   P_mat(i,j) * Conv_lim];
                Constraints = [Constraints, I_C(i, j, k) >= - P_mat(i,j) * Conv_lim];
                Constraints = [Constraints, I_C(i, j, k) == - I_C(j, i, k)];
                end
            end
        end
    else
        for i = 1:Node_num
            for j = i:Node_num
                if(P_mat(i,j)==1)
                    for k = 1:T_num
                        Constraints = [Constraints, I_C(i,j,k) == - I_C(j,i,k)];
                    end
                else
                    for k = 1:T_num
                        Constraints = [Constraints, I_C(i,j,k)==0];
                        Constraints = [Constraints, I_C(j,i,k)==0];
                    end
                end
            end
        end

        for k = 1:T_num
            Constraints = [Constraints,sum(sum(abs(I_C(:,:,k))))/2 <= P_lim * lim_sum];
        end
    end

    % [4] the constriant for battery Power(current) Limit
    for i = 1:Node_num
        for k = 1:T_num
        Constraints = [Constraints, SYS.Bat{i}.curlim >= I_B(i, k)];
        Constraints = [Constraints, I_B(i, k) >= -SYS.Bat{i}.curlim];
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
                Constraints = [Constraints, I_C(i, j, k) == I_C(i, j, 1)];
                end
            end
        end
    end

    % [7] the constriant for output voltage and current

    for k = 1:T_num
        Constraints = [Constraints, I_L(k) == (abs(floor(Output_waveform(k) /24))/Node_num)*I_L(Node_num+1)];
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
   
    % get the optimal value for the first optimization
    Optimal_IL = value(sum(I_L));

    yalmip('clear');
    % Second optimization - system efficiency
    ops = sdpsettings('solver', 'gurobi', ...    
                      'gurobi.Threads', 8, ...
                      'gurobi.MIPGap', 0.01, ... 
                      'gurobi.TimeLimit', 120, ...
                      'gurobi.Presolve', 2, ...
                      'gurobi.Cuts', 2,...
                      'gurobi.Heuristics', 1);

    I_B = sdpvar(Node_num, T_num, 'full');
    I_L = sdpvar(T_num, 1, 'full');
    Q_B = sdpvar(Node_num, T_num, 'full');
    Q_C = sdpvar(Node_num, T_num, 'full');
    I_C = sdpvar(Node_num, Node_num, T_num, 'full');
    Bat_use = QL_connection;
  
    if (TEST.Conv.flow == 0)
    I_process = 0;
    for k = 1:T_num
        I_process = I_process + sum(sum(abs(I_C(:,:,k))))/2 * Delta_t(k);
    end
    else
         I_process = sum(sum(abs(I_C(:,:,1))))/2;
    end

    Objective2 = I_process; % minimize converter processing power

    Constraints = [];
    
    %% constraints
    % maintains maximum output
    Constraints = [Constraints, sum(I_L) == Optimal_IL];

    % [1] the constraint for charge matrix
    for i = 1:Node_num
        for j = 1:T_num
        Constraints = [Constraints, Q_B(i, j) == I_B(i, j) * Delta_t(j)];
        Constraints = [Constraints, Q_C(i, j) == sum(I_C(i, :, j)) * Delta_t(j)];
        end
    end

    % [2] the constriant for KCL 
    M = 100;
    for k = 1:T_num
         for i = 1:Node_num
         Constraints = [Constraints, sum(I_C(i, :, k)) + I_B(i,k) - I_L(k) <= M * (1 - Bat_use(i, k))];
         Constraints = [Constraints, I_L(k) - sum(I_C(i, :, k)) - I_B(i,k) <= M * (1 - Bat_use(i, k))];
         end
    end

    % [3] the constriant for differential power processing links and Converter Power(current) Limit

    if(TEST.Conv.partition==1)
        for i = 1:Node_num
            for j = i:Node_num
                for k = 1:T_num
                Constraints = [Constraints, I_C(i, j, k) <=   P_mat(i,j) * Conv_lim];
                Constraints = [Constraints, I_C(i, j, k) >= - P_mat(i,j) * Conv_lim];
                Constraints = [Constraints, I_C(i, j, k) == - I_C(j, i, k)];
                end
            end
        end
    else
        for i = 1:Node_num
            for j = i:Node_num
                if(P_mat(i,j)==1)
                    for k = 1:T_num
                        Constraints = [Constraints, I_C(i,j,k) == - I_C(j,i,k)];
                    end
                else
                    for k = 1:T_num
                        Constraints = [Constraints, I_C(i,j,k)==0];
                        Constraints = [Constraints, I_C(j,i,k)==0];
                    end
                end
            end
        end

        for k = 1:T_num
            Constraints = [Constraints,sum(sum(abs(I_C(:,:,k))))/2 <= P_lim * lim_sum];
        end
    end

    % [4] the constriant for battery Power(current) Limit
    for i = 1:Node_num
        for k = 1:T_num
        Constraints = [Constraints, SYS.Bat{i}.curlim >= I_B(i, k)];
        Constraints = [Constraints, I_B(i, k) >= -SYS.Bat{i}.curlim];
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
                Constraints = [Constraints, I_C(i, j, k) == I_C(i, j, 1)];
                end
            end
        end
    end

    % [7] the constriant for output voltage and current
    for k = 1:T_num
        Constraints = [Constraints, I_L(k) == (abs(floor(Output_waveform(k) /24))/Node_num)*I_L(Node_num+1)];
    end
        % battery connect to the load
        M = 100;
        for k = 1:T_num
            for i = 1:Node_num
            Constraints = [Constraints, I_B(i, k) <= M * Bat_use(i, k)];
            Constraints = [Constraints, I_B(i, k) >= -M * Bat_use(i, k)];
            end
        end

    optimize(Constraints, Objective2, ops);

    OptRes.I_process = value(I_process);
    OptRes.I_B = value(I_B);
    OptRes.I_L = value(I_L);
    OptRes.I_C = value(I_C);
    OptRes.Q_B = value(Q_B);
    OptRes.Q_C = value(Q_C);
    OptRes.SOC = value(SOC);
    OptRes.Bat_use = value(Bat_use);
    OptRes.up = value(I_L(Node_num+1))*Node_num/lim_sum;

    % calculate efficiency
    for i = 1:T_num
        energy_t(i) = Delta_t(i) * value(I_L(i));
    end
    energy_sum = sum(energy_t)*Node_num;
    conv_eta = 0.85;
    OptRes.efficience = (1-(OptRes.I_process*(1-conv_eta))/energy_sum)*100;

end
