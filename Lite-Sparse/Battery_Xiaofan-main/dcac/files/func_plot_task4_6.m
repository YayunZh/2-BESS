function y = func_plot_task4_6(ENV,ExpTrajRes,ExpTrajRes_trad)

for ind = 1:size(ExpTrajRes,4)
    u_gap(ind) = max(ExpTrajRes{1,1,1,ind}{1,1}.uratio_energy.value) - max(ExpTrajRes_trad{1,1,1,ind}{1,1}.uratio_energy.value);
end
[u_gap_max, u_gap_ind] = max(u_gap);
ExpTrajRes_typ = ExpTrajRes{1,1,1,u_gap_ind}{1,1};
ExpTrajRes_trad_typ = ExpTrajRes_trad{1,1,u_gap_ind}{1,1};

% ExpTrajRes_typ.Ebessout = ExpectedTrajEbessout (ExpTrajRes, ENV);
% ExpTrajRes_trad_typ.Ebessout = ExpectedTrajEbessout (ExpTrajRes_trad, ENV);

scalefac1 = 150/36;

%% Typical Trajectory Plots
scalefac1 = 150/36;
fig1 = figure();
left_color = [0.8500, 0.3250, 0.0980];
right_color = [0 0.4470 0.7410];
set(fig1,'defaultAxesColorOrder',[left_color; right_color]);
%colororder({'b','m'});

yyaxis right
l_pbess_out = stairs(ExpTrajRes_typ.Pbessout.abs_time(1:end-1), scalefac1*ExpTrajRes_typ.Pbessout.value(2:end),':','color',[0.4940, 0.1840, 0.5560],'linewidth',1.5);
hold on;
l_pgrid = plot(0:0.1:24, -scalefac1*ENV.Pgrid_fit(0:0.1:24),':b','linewidth',1.5);
hold on;
l_pevin = stairs(ExpTrajRes_typ.Pevin.abs_time(1:end-1), scalefac1*ExpTrajRes_typ.Pevin.value(2:end),':','color',[0.4660, 0.6740, 0.1880],'linewidth',1.5);
%hold on;
%l_pgrid_avg = stairs(ExpTrajRes_typ.Pgrid_avg.abs_time(1:end-1), -scalefac1*ExpTrajRes_typ.Pgrid_avg.value(2:end),'color',[0.4660, 0.6740, 0.1880],'linewidth',1.5);
ylim([-60,160]);
ylabel('Power (kW)');
hold on;
% l_ev_demand = stairs(ExpTrajRes_typ.Evdemand.abs_time(1:end-1), scalefac1*ExpTrajRes_typ.Evdemand.value(2:end),'linewidth',1.5);
% hold on;
% l_ev_demand = stairs(ExpTrajRes_typ.Tdelay.abs_time(1:end-1), scalefac1*ExpTrajRes_typ.Tdelay.value(2:end),'linewidth',1.5);
yyaxis left
l_ebess = plot(ExpTrajRes_typ.Ebess.abs_time, scalefac1*ExpTrajRes_typ.Ebess.value,'color',[0.8500, 0.3250, 0.0980],'linewidth',2.5);
ylim([-30,80]);
ylabel('Energy (kWh)');
xlim([0,25]);
xlabel('Time (h)');
title('LS-HiPPP');
%legend([l_pgrid_avg, l_pgrid, l_pbess_out, l_ebess, l_pevin],'Avg Power from Grid (kW)','Grid Constraint (kW)','BESS Power (kW)','BESS Energy (kWh)','EV Charging Power (kW)');
legend([l_pgrid, l_pbess_out, l_ebess, l_pevin],'Grid Power Constraint','BESS Power','BESS Energy','EV Charging Power');

fig2 = figure();
set(fig2,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis right
l_pbess_out_trad = stairs(ExpTrajRes_trad_typ.Pbessout.abs_time(1:end-1), scalefac1*ExpTrajRes_trad_typ.Pbessout.value(2:end),':','color',[0.4940, 0.1840, 0.5560],'linewidth',1.5);
hold on;
l_pgrid_trad = plot(0:0.1:24, -scalefac1*ENV.Pgrid_fit(0:0.1:24),':b','linewidth',1.5);
hold on;
l_pevin_trad = stairs(ExpTrajRes_trad_typ.Pevin.abs_time(1:end-1), scalefac1*ExpTrajRes_trad_typ.Pevin.value(2:end),':','color',[0.4660, 0.6740, 0.1880],'linewidth',1.5);
%hold on;
%l_pgrid_avg = stairs(ExpTrajRes_trad_typ.Pgrid_avg.abs_time(1:end-1), -scalefac1*ExpTrajRes_trad_typ.Pgrid_avg.value(2:end),'color',[0.4660, 0.6740, 0.1880],'linewidth',1.5);
ylim([-60,160]);
ylabel('Power (kW)');
hold on;
% l_ev_demand = stairs(ExpTrajRes_trad_typ.Evdemand.abs_time(1:end-1), scalefac1*ExpTrajRes_trad_typ.Evdemand.value(2:end),'linewidth',1.5);
% hold on;
% l_ev_demand = stairs(ExpTrajRes_trad_typ.Tdelay.abs_time(1:end-1), scalefac1*ExpTrajRes_trad_typ.Tdelay.value(2:end),'linewidth',1.5);
yyaxis left
l_ebess_trad = plot(ExpTrajRes_trad_typ.Ebess.abs_time, scalefac1*ExpTrajRes_trad_typ.Ebess.value,'color',[0.8500, 0.3250, 0.0980],'linewidth',2.5);
ylim([-30,80]);
ylabel('Energy (kWh)');
xlim([0,25]);
xlabel('Time (h)');

title('C-PPP');
%legend([l_pgrid_avg, l_pgrid, l_pbess_out, l_ebess, l_pevin],'Avg Power from Grid (kW)','Grid Constraint (kW)','BESS Power (kW)','BESS Energy (kWh)','EV Charging Power (kW)');
legend([l_pgrid_trad, l_pbess_out_trad, l_ebess_trad, l_pevin_trad],'Grid Power Constraint','BESS Power','BESS Energy','EV Charging Power');

%% Grid Constraint and EV Demand Gaps
fig3 = figure();
set(fig3,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis right
l_pgrid = plot(0:0.1:24, scalefac1*ENV.Pgrid_fit(0:0.1:24),'linewidth',1.5,'color','b');
ylabel('Power (kW)');
ylim([15,60]);
grid on;
grid minor;
hold on;

yyaxis left
l_ev_demand = stairs(ExpTrajRes_typ.Evdemand.abs_time(1:end-1), scalefac1*ExpTrajRes_typ.Evdemand.value(2:end),'linewidth',1.5);
hold on;
%l_grid_energy = stairs(ExpTrajRes_typ.Evdemand.abs_time(1:end-1), scalefac1*ExpTrajRes_typ.L3ChargeTime.value(2:end).*ExpTrajRes_typ.Pgrid_avg.value(2:2:end),'color',[0.4660, 0.6740, 0.1880],'linewidth',1.5);
l_grid_energy = stairs(ExpTrajRes_typ.Evdemand.abs_time(1:end-1), scalefac1*ExpTrajRes_typ.L3ChargeTime.value(2:end).*ExpTrajRes_typ.Pgrid_avg.value(2:2:end),':','linewidth',1.5);
ylabel('Energy (kWh)');
%legend([l_pgrid,l_ev_demand,l_grid_energy],'Grid Power Constraint','EV Energy Demand','Grid Energy');
xlim([0,25]);
ylim([0,50]);
xlabel('Time (h)');
grid on;
grid minor;

%% Measure of Dispersions Utilization
% figure();
% [ExpTrajRes_avg.uratio_energy,~] = ExpectedTrajuratio_energy(ExpTrajRes, ENV);
% [ExpTrajRes_trad_avg.uratio_energy,~] = ExpectedTrajuratio_energy(ExpTrajRes_trad, ENV);
% l_uratio_energy = plot(ExpTrajRes_avg.uratio_energy.abs_time(1:end-1), 100*ExpTrajRes_avg.uratio_energy.prcgap(2:end),...
%     '-s','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue','color','b','linewidth',2);
% hold on;
% l_uratio_energy_trad = plot(ExpTrajRes_trad_avg.uratio_energy.abs_time(1:end-1), 100*ExpTrajRes_trad_avg.uratio_energy.prcgap(2:end),...
%     '-d','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','color','r','linewidth',2);
% xlabel('Time (h)');
% ylabel('Interdecitile Range of Battery Capacity Utilization (%)')
% legend([l_uratio_energy,l_uratio_energy_trad],'LS-HiPPP','PPP');
% grid on;
% grid minor;
% 
% figure();
% [ExpTrajRes_avg.uratio_energy,~] = ExpectedTrajuratio_energy(ExpTrajRes, ENV);
% [ExpTrajRes_trad_avg.uratio_energy,~] = ExpectedTrajuratio_energy(ExpTrajRes_trad, ENV);
% l_uratio_energy = plot(ExpTrajRes_avg.uratio_energy.abs_time(1:end-1), 100*3*ExpTrajRes_avg.uratio_energy.std(2:end),...
%     '-s','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue','color','b','linewidth',2);
% hold on;
% l_uratio_energy_trad = plot(ExpTrajRes_trad_avg.uratio_energy.abs_time(1:end-1), 100*3*ExpTrajRes_trad_avg.uratio_energy.std(2:end),...
%     '-d','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','color','r','linewidth',2);
% xlabel('Time (h)');
% ylabel('3\sigma Performance Spread of Battery Capacity Utilization (%)')
% legend([l_uratio_energy,l_uratio_energy_trad],'LS-HiPPP','PPP');
% grid on;
% grid minor;

%% Measure of Dispersions Energy Out (Normalized)
figure();
[ExpTrajRes_avg.Ebessout,~] = ExpectedTrajEbessout(ExpTrajRes, ENV);
[ExpTrajRes_trad_avg.Ebessout,~] = ExpectedTrajEbessout(ExpTrajRes_trad, ENV);
l_Ebessout = plot(ExpTrajRes_avg.Ebessout.abs_time(1:end-1), 100*ExpTrajRes_avg.Ebessout.prcgap(2:end)./ExpTrajRes_avg.Ebessout.mean(2:end),...
    '-s','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue','color','b','linewidth',2);
hold on;
l_Ebessout_trad = plot(ExpTrajRes_trad_avg.Ebessout.abs_time(1:end-1), 100*ExpTrajRes_trad_avg.Ebessout.prcgap(2:end)./ExpTrajRes_avg.Ebessout.mean(2:end),...
    '-d','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','color','r','linewidth',2);
xlabel('Time (h)');
ylabel('IDR of 2-BESS Energy Output (%)');
legend([l_Ebessout,l_Ebessout_trad],'LS-HiPPP','C-PPP');
grid on;
grid minor;

figure();
[ExpTrajRes_avg.Ebessout,~] = ExpectedTrajEbessout(ExpTrajRes, ENV);
[ExpTrajRes_trad_avg.Ebessout,~] = ExpectedTrajEbessout(ExpTrajRes_trad, ENV);
l_Ebessout = plot(ExpTrajRes_avg.Ebessout.abs_time(1:end-1), 100*3*ExpTrajRes_avg.Ebessout.std(2:end)./ExpTrajRes_avg.Ebessout.mean(2:end),...
    '-s','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue','color','b','linewidth',2);
hold on;
l_Ebessout_trad = plot(ExpTrajRes_trad_avg.Ebessout.abs_time(1:end-1), 100*3*ExpTrajRes_trad_avg.Ebessout.std(2:end)./ExpTrajRes_avg.Ebessout.mean(2:end),...
    '-d','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','color','r','linewidth',2);
xlabel('Time (h)');
ylabel('3\sigma spread of 2-BESS Energy Output (%)');
legend([l_Ebessout,l_Ebessout_trad],'LS-HiPPP','C-PPP');
grid on;
grid minor;

%% Measure of Dispersions Energy Out (True Values)
% figure();
% [ExpTrajRes_avg.Ebessout,~] = ExpectedTrajEbessout(ExpTrajRes, ENV);
% [ExpTrajRes_trad_avg.Ebessout,~] = ExpectedTrajEbessout(ExpTrajRes_trad, ENV);
% l_Ebessout = plot(ExpTrajRes_avg.Ebessout.abs_time(1:end-1), scalefac1*ExpTrajRes_avg.Ebessout.prcgap(2:end),...
%     '-s','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue','color','b','linewidth',2);
% hold on;
% l_Ebessout_trad = plot(ExpTrajRes_trad_avg.Ebessout.abs_time(1:end-1), scalefac1*ExpTrajRes_trad_avg.Ebessout.prcgap(2:end),...
%     '-d','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','color','r','linewidth',2);
% xlabel('Time (h)');
% ylabel('IDR of 2-BESS Energy Output (%)');
% legend([l_Ebessout,l_Ebessout_trad],'LS-HiPPP','C-PPP');
% grid on;
% grid minor;
% 
% figure();
% [ExpTrajRes_avg.Ebessout,~] = ExpectedTrajEbessout(ExpTrajRes, ENV);
% [ExpTrajRes_trad_avg.Ebessout,~] = ExpectedTrajEbessout(ExpTrajRes_trad, ENV);
% l_Ebessout = plot(ExpTrajRes_avg.Ebessout.abs_time(1:end-1), scalefac1*3*ExpTrajRes_avg.Ebessout.std(2:end),...
%     '-s','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue','color','b','linewidth',2);
% hold on;
% l_Ebessout_trad = plot(ExpTrajRes_trad_avg.Ebessout.abs_time(1:end-1), scalefac1*3*ExpTrajRes_trad_avg.Ebessout.std(2:end),...
%     '-d','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','color','r','linewidth',2);
% xlabel('Time (h)');
% ylabel('3\sigma spread of 2-BESS Energy Output (%)');
% legend([l_Ebessout,l_Ebessout_trad],'LS-HiPPP','C-PPP');
% grid on;
% grid minor;

%% Utilization Boxplot
[ExpTrajRes_typ.uratio_energy,uratio_energy_matrix] = ExpectedTrajuratio_energy(ExpTrajRes, ENV);
[ExpTrajRes_trad_typ.uratio_energy,uratio_energy_trad_matrix] = ExpectedTrajuratio_energy(ExpTrajRes_trad, ENV);
figure();
l_uratio_energy2 = boxplot(100*uratio_energy_matrix, ...
    'positions',round(ExpTrajRes_typ.uratio_energy.abs_time,1), ...
    'labels', round(ExpTrajRes_typ.uratio_energy.abs_time,1));
xtickangle(45);
xlabel('Time (h)');
ylabel('Battery Energy Utilization (%)');
ylim([50,100]);
title('LS-HiPPP');
% box_plot_1 = findobj(gca,'Tag','Box');
% legend([box_plot_1(13), box_plot_1(1)],{'LS-HiPPP','Conventional PPP'});
set(findobj(gca,'type','line'),'linew',1.5);
saveas(gcf,'Tasl4_6_1.fig');

figure();
l_uratio_energy2_trad = boxplot(100*uratio_energy_trad_matrix,...
    'positions', round(ExpTrajRes_trad_typ.uratio_energy.abs_time,1), ....
    'labels', round(ExpTrajRes_trad_typ.uratio_energy.abs_time,1));
xtickangle(45);
xlabel('Time (h)');
ylabel('Battery Energy Utilization (%)');
ylim([50,100]);
title('C-PPP');
% box_plot_1 = findobj(gca,'Tag','Box');
% legend([box_plot_1(13), box_plot_1(1)],{'LS-HiPPP','Conventional PPP'});
set(findobj(gca,'type','line'),'linew',1.5);
saveas(gcf,'Tasl4_6_2.fig');

%% User Defined Functions
function [ExpTrajRes_uratio_energy, uratio_energy_matrix] = ExpectedTrajuratio_energy(ExpTrajRes, ENV)
    
    uratio_energy_matrix = [];
    for i = 1:ENV.Var_Conv.MC_trial
        uratio_energy_matrix = [uratio_energy_matrix; ExpTrajRes{1,1,i}{1,1}.uratio_energy.value];
    end

    ExpTrajRes_uratio_energy.delta_time = ExpTrajRes{1,1,1}{1,1}.uratio_energy.delta_time;
    ExpTrajRes_uratio_energy.abs_time = ExpTrajRes{1,1,1}{1,1}.uratio_energy.abs_time;
    ExpTrajRes_uratio_energy.upper = max(uratio_energy_matrix,[],1);
    ExpTrajRes_uratio_energy.lower = min(uratio_energy_matrix,[],1);
    ExpTrajRes_uratio_energy.mean = mean(uratio_energy_matrix,1);
    ExpTrajRes_uratio_energy.std = std(uratio_energy_matrix);
    ExpTrajRes_uratio_energy.prcgap = prctile(uratio_energy_matrix,90) - prctile(uratio_energy_matrix,10);

end

function [ExpTrajRes_Ebessout,Ebessout_matrix] = ExpectedTrajEbessout(ExpTrajRes, ENV)
    
    Ebessout_matrix = [];
    for i = 1:ENV.Var_Conv.MC_trial
        Ebessout_matrix = [Ebessout_matrix; ExpTrajRes{1,1,i}{1,1}.Ebessout.value];
    end
    ExpTrajRes_Ebessout.delta_time = ExpTrajRes{1,1,1}{1,1}.uratio_energy.delta_time;
    ExpTrajRes_Ebessout.abs_time = ExpTrajRes{1,1,1}{1,1}.uratio_energy.abs_time;
    ExpTrajRes_Ebessout.upper = max(Ebessout_matrix,[],1);
    ExpTrajRes_Ebessout.lower = min(Ebessout_matrix,[],1);
    ExpTrajRes_Ebessout.mean = mean(Ebessout_matrix,1);
    ExpTrajRes_Ebessout.std = std(Ebessout_matrix);
    ExpTrajRes_Ebessout.prcgap = prctile(Ebessout_matrix,90) - prctile(Ebessout_matrix,10);
    

end


end

