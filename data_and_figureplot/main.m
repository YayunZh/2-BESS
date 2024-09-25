% 
% %% Part 1 comparision of dc-ac and mac-2BESS
% % figure 1 mac-2BESS boxplot
%     load("MC5%.mat")
%     heter = 1; % heterogeniety 1-5% 2-10% 3-15% 4-20%
%     num_converterrating = 18;
%     num_MC = 10;
% 
%     for i=1:num_converterrating
%         for j=1:num_MC
%             MCEff(i+4,j) = Opt_MC.MC(heter,i+1,j).efficience2;
%             MCUp(i+4,j) = 100*Opt_MC.MC(heter,i+1,j).up;
%         end
%     end
% 
%     load("MC0.1.mat")
%     heter = 1; % heterogeniety 1-5% 2-10% 3-15% 4-20%
%     num_converterrating = 4;
%     num_MC = 10;
% 
%     for i=1:num_converterrating
%         for j=1:num_MC
%             MCEff(i,j) = Opt_MC.MC(heter,i,j).efficience2;
%             MCUp(i,j) = 100*Opt_MC.MC(heter,i,j).up;
%         end
%     end
% 
% 
%     % battery power utilization
%     func_figure1(heter,MCUp);
%     % system power efficiency
%     func_figure1_2(heter,MCEff);


% % figure 2 battery power utilization comparision
%     load("MC20%.mat");
%     heter = 4; % heterogeniety 1-5% 2-10% 3-15% 4-20%
%     num_converterrating = 18;
% 
%     for i = 1:num_converterrating
%     mac(1,i) = 100*Opt_MC.power_Utilization(heter,i);
%     end
%     func_figure2(heter,mac);
% 
% % figure 3 system power efficiency comparision
%     load("MC20%.mat");
%     heter = 4; % heterogeniety 1-5% 2-10% 3-15% 4-20%
%     num_converterrating = 18;
% 
%     for i = 1:num_converterrating
%     mac(1,i) = Opt_MC.efficience2(heter,i);
%     end
%     func_figure3(heter,mac);

% Part 2 PQcurve
% figure 4 mac-2BESS PQ
% load("PQI_max10%.mat");
% heter = 2;
% for i=1:10
%     for j=1:100
%         Ipeak(i,j)=Power_capability(heter,i,j)*0.98;
%     end
% end
% func_figure4(heter,Ipeak);

% figure 5 DC-AC PQ
% load("dcacEff.mat");
% load("dcacUp.mat");
% heter = 2;
% for i=1:10
%     Ipeak(1,i)= 0.0025*dcacEff(heter,i)*dcacUp(heter,i);
% end
% func_figure5(heter,Ipeak);

% figure 6 fitting line
% 定义多边形的顶点
x = [2100, 5100, 2100, -2100,-5300, -2100,2100];
y = [2700, 0, -2700, -2700, 0,2700,2700];

func_figure6(x,y)














