% Sb=100MVA,Vb=12.66KV
clear 
clc
tic
warning off
load ("network_data.mat")

%% SYSTEM Parameter
mpc = IEEE33BW; % IEEE 33-bus system 
pload = (data(:,5))/1e5; % Pload;
qload = (data(:,6))/1e5; % Qload;
branch = mpc.branch;
branch(:,3) = branch(:,3)*100/(12.66^2); % Calculate per-unit impedance
r = real(branch(:,3));
x = imag(branch(:,3));
T = 1;
nb = 33; 
nl = 32; 

upstream = zeros(nb,nl);
dnstream = zeros(nb,nl);
for i = 1:nl
    upstream(i,i) = 1;
end
for i=[1:16,18:20,22:23,25:31]
    dnstream(i,i+1) = 1;
end
dnstream(1,18) = 1;
dnstream(2,22) = 1;
dnstream(5,25) = 1;
dnstream(33,1) = 1;

% voltage profile
Vmax = [1*1*ones(nb-1,T)
        1*1*ones(1,T)];
Vmin = [0.93*0.93*ones(nb-1,T)
        1*1*ones(1,T)];
Pgmax = [zeros(nb-1,T)
         ones(1,T)];
Qgmax = [zeros(nb-1,T)
         ones(1,T)];

%% The number of buses for installing PVs or batteries

PVnum = 2;
BESSnum = 2;

%% Part 1 DC-AC
%% variables
V = sdpvar(nb,T); % Voltage squared 
I = sdpvar(nl,T); % current squared
P = sdpvar(nl,T); % line active power
Q = sdpvar(nl,T); % line reactive power
Pg = sdpvar(nb,T); % generator active power
Qg = sdpvar(nb,T); % generator reactive power

PVsite = binvar(nb, T, 'full');
Ppv = sdpvar(nb, T, 'full');
Qpv = sdpvar(nb, T, 'full');

BESSsite = binvar(nb, T, 'full');
Pbess = sdpvar(nb, T, 'full');
Qbess = sdpvar(nb, T, 'full');

%% constraints
Constraints = [];

% Bus power constraints
Pin = -upstream * P + upstream * (I.*(r*ones(1,T))) + dnstream * P - BESSsite .* Pbess - PVsite .* Ppv ;   % injected active power
Qin = -upstream * Q + upstream * (I.*(x*ones(1,T))) + dnstream * Q - BESSsite .* Qbess - PVsite .* Qpv ;   % injected reactive power
Constraints = [Constraints, Pin + pload - Pg == 0];
Constraints = [Constraints, Qin + qload - Qg == 0];

% Distributed energy constraints
% PV
Constraints = [Constraints,  0 <= Ppv, sum(Ppv(:,1)) == 0.01];
Constraints = [Constraints,  0 == Qpv];
Constraints = [Constraints,  0 <= Pbess];

Constraints = [Constraints, sum(PVsite) == PVnum];
Constraints = [Constraints, PVsite(33) == 0];
Constraints = [Constraints, PVsite(1) == 0];

Constraints = [Constraints, sum(BESSsite) == BESSnum];
Constraints = [Constraints, BESSsite(33) == 0];
Constraints = [Constraints, BESSsite(1) == 0];

M=0.01;
for i = 1:nb
    
    Constraints = [Constraints,  0.005 >= Pbess(i)];
    Constraints = [Constraints, Pbess(i)^2 + Qbess(i)^2 <= (0.005*0.005)];

    Constraints = [Constraints, Qpv(i) <= M * PVsite(i), - M * PVsite(i) <= Qpv(i)];
    Constraints = [Constraints, Ppv(i) <= M * PVsite(i), - M * PVsite(i) <= Ppv(i)];

    Constraints = [Constraints, Qbess(i) <= M * BESSsite(i), - M * BESSsite(i) <= Qbess(i)];
    Constraints = [Constraints, Pbess(i) <= M * BESSsite(i), - M * BESSsite(i) <= Pbess(i)];

end

% 
Constraints = [Constraints, V(branch(:,2),:) == V(branch(:,1),:) - 2*(r*ones(1,T)).*P - 2*(x*ones(1,T)).*Q + ((r.^2+x.^2)*ones(1,T)).*I];

% Second-order cone constraints
 for t=1:T
 for i = 1:32
     Constraints = [Constraints, cone([2*P(i,t);2*Q(i,t);I(i,t)-V(branch(i, 1),t)],I(i,t)+V(branch(i, 1),t))];
 end
 end

Constraints = [Constraints, Vmin <= V, V <= Vmax];
Constraints = [Constraints, -Pgmax <= Pg, Pg <= Pgmax, -Qgmax <= Qg, Qg <= Qgmax];
Constraints = [Constraints, -0.21 <= I, I <= 0.21];
Constraints = [Constraints, -0.21 <= P, P <= 0.11, -0.21 <= Q, Q <= 0.21];

%% objective function
objective = sum(sum(I.*(r*ones(1,T)))) + sum(sum(I.*(x*ones(1,T))));
toc 

%%
ops=sdpsettings('verbose', 1, 'solver', 'gurobi');
sol=optimize(Constraints,objective,ops);
toc 
lineloss = 100*(sum(sum(I.*(r*ones(1,T)))) + sum(sum(I.*(x*ones(1,T)))));
%%
if sol.problem == 0
    disp('succcessful solved');
else
    disp('error');
    yalmiperror(sol.problem)
end
%%
 ResOpt.V = sqrt(value(V));
 ResOpt.I = value(I);
 ResOpt.P = value(P);
 ResOpt.Q = value(Q);
 ResOpt.Pg = value(Pg);
 ResOpt.Qg = value(Qg);
 ResOpt.Wtsite = value(PVsite);
 ResOpt.Pw = value(Ppv);
 ResOpt.Qw = value(Qpv);
 ResOpt.BESSsite = value(BESSsite);
 ResOpt.Pbess = value(Pbess);
 ResOpt.Qbess = value(Qbess);
 ResOpt.Loss = value(lineloss);
 
 % save result
 datapro(1,1) = 1;
    for b = 1:32
        datapro(1,b+1) = ResOpt.V(b,1);
    end
    loss(1) = (1-value(lineloss)/0.1844)*100;
    V_avg(1) = sum(ResOpt.V)/33;

 yalmip('clear');

%% Part 2 MAC-2BESS

%% variables
V = sdpvar(nb,T); % Voltage squared 
I = sdpvar(nl,T); % current squared
P = sdpvar(nl,T); % line active power
Q = sdpvar(nl,T); % line reactive power
Pg = sdpvar(nb,T); % generator active power
Qg = sdpvar(nb,T); % generator reactive power

PVsite = binvar(nb, T, 'full');
Ppv = sdpvar(nb, T, 'full');
Qpv = sdpvar(nb, T, 'full');

BESSsite = binvar(nb, T, 'full');
Pbess = sdpvar(nb, T, 'full');
Qbess = sdpvar(nb, T, 'full');

%% constraints
Constraints = [];

Pin = -upstream * P + upstream * (I.*(r*ones(1,T))) + dnstream * P - BESSsite .* Pbess - PVsite .* Ppv ;   
Qin = -upstream * Q + upstream * (I.*(x*ones(1,T))) + dnstream * Q - BESSsite .* Qbess - PVsite .* Qpv ;   
Constraints = [Constraints, Pin + pload - Pg == 0];
Constraints = [Constraints, Qin + qload - Qg == 0];

Constraints = [Constraints,  0 <= Ppv, sum(Ppv(:,1)) == 0.01];
Constraints = [Constraints,  0 == Qpv];
Constraints = [Constraints,  0 <= Pbess];

Constraints = [Constraints, sum(PVsite) == PVnum];
Constraints = [Constraints, PVsite(33) == 0];
Constraints = [Constraints, PVsite(1) == 0];

Constraints = [Constraints, sum(BESSsite) == BESSnum];
Constraints = [Constraints, BESSsite(33) == 0];
Constraints = [Constraints, BESSsite(1) == 0];

M=0.01;
for i = 1:nb
    
    Constraints = [Constraints,  0.005*1.89 >= Qbess(i)];
    Constraints = [Constraints,  0.005 >= Pbess(i)];

    Constraints = [Constraints, Pbess(i) <= (-0.9*(Qbess(i)/0.005)+1.7)*0.005];
    Constraints = [Constraints, Qpv(i) <= M * PVsite(i), - M * PVsite(i) <= Qpv(i)];
    Constraints = [Constraints, Ppv(i) <= M * PVsite(i), - M * PVsite(i) <= Ppv(i)];

    Constraints = [Constraints, Qbess(i) <= M * BESSsite(i), - M * BESSsite(i) <= Qbess(i)];
    Constraints = [Constraints, Pbess(i) <= M * BESSsite(i), - M * BESSsite(i) <= Pbess(i)];

end

Constraints = [Constraints, V(branch(:,2),:) == V(branch(:,1),:) - 2*(r*ones(1,T)).*P - 2*(x*ones(1,T)).*Q + ((r.^2+x.^2)*ones(1,T)).*I];

 for t=1:T
 for i = 1:32
     Constraints = [Constraints, cone([2*P(i,t);2*Q(i,t);I(i,t)-V(branch(i, 1),t)],I(i,t)+V(branch(i, 1),t))];
 end
 end

Constraints = [Constraints, Vmin <= V, V <= Vmax];
Constraints = [Constraints, -Pgmax <= Pg, Pg <= Pgmax, -Qgmax <= Qg, Qg <= Qgmax];
Constraints = [Constraints, -0.21 <= I, I <= 0.21];
Constraints = [Constraints, -0.21 <= P, P <= 0.11, -0.21 <= Q, Q <= 0.21];

%% objective function
objective = sum(sum(I.*(r*ones(1,T)))) + sum(sum(I.*(x*ones(1,T))));

toc 

ops=sdpsettings('verbose', 1, 'solver', 'gurobi');
sol=optimize(Constraints,objective,ops);
toc 
lineloss = 100*(sum(sum(I.*(r*ones(1,T)))) + sum(sum(I.*(x*ones(1,T)))));

%% 
if sol.problem == 0
    disp('succcessful solved');
else
    disp('error');
    yalmiperror(sol.problem)
end

 ResOpt.V = sqrt(value(V));
 ResOpt.I = value(I);
 ResOpt.P = value(P);
 ResOpt.Q = value(Q);
 ResOpt.Pg = value(Pg);
 ResOpt.Qg = value(Qg);
 ResOpt.Wtsite = value(PVsite);
 ResOpt.Pw = value(Ppv);
 ResOpt.Qw = value(Qpv);
 ResOpt.BESSsite = value(BESSsite);
 ResOpt.Pbess = value(Pbess);
 ResOpt.Qbess = value(Qbess);
 ResOpt.Loss = value(lineloss);
 
 % save results
 datapro(2,1) = 1;
    for b = 1:32
        datapro(2,b+1) = ResOpt.V(b,1);
    end

    loss(2) = (1-value(lineloss)/0.1844)*100;
    V_avg(2) = sum(ResOpt.V)/33;
   
x = linspace(1,33,33);

save("power loss minimization","loss");
save("average voltage","V_avg");


%% figure plot
figure;
hold on;
set(gcf, 'Color', 'w');

num_colors = 9;
linecolors = [linspace(0.9, 0.1, num_colors)', linspace(0.8, 0.2, num_colors)',  linspace(1,0.7, num_colors)'];

plot(x, datapro(1,:), '-o', 'MarkerFaceColor', [0.4	0.8	1], 'MarkerSize', 1, 'LineWidth', 1, 'DisplayName', 'DC-AC');
plot(x, datapro(2,:), '-o', 'MarkerFaceColor', linecolors(2, :), 'MarkerSize', 1, 'LineWidth', 1, 'DisplayName', 'MAC-2BESS');

legend('show');

xlim([1 35]);
ylim([0.95 1]); 

title('Buses voltage profile');
xlabel('Bus Number');
ylabel('Voltage');

grid on;

hold off;
