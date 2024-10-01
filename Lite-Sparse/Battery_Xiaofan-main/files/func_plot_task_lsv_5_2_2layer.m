function y = func_plot_task_lsv_5_2_2layer(CL,ENV,ExpRes,ENV_trad,ExpRes_trad,ENV_fpp,ExpRes_fpp)
    %% Plot the Utilization -- Initial design Converters VS Traditional
    figure(15);
        for i = 1:size(ExpRes.bat_uratio_energy,1)
            subplot(4,2,i);
            pic(i) = plot((ENV.Fir_Conv.Num*ENV.Sweep.Conv.e_lim)/(CL.Stat.Bat_num*CL.Bat{1}.curlim*CL.Bat{1}.volt),...
            100*mean(reshape(ExpRes.bat_uratio_energy(i,:,:),ENV.Sweep.Stat.Conv,ENV.Fir_Conv.MC_trial),2),'d-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_energy,1),0,0]);
            hold on;
            pic_trad(i) = plot(ENV_trad.Sweep.Conv_energy_sum/ENV_trad.Sweep.Bat_energy_sum, ...
            100*mean(reshape(ExpRes_trad.bat_uratio_energy(i,:,:),ENV_trad.Sweep.Stat.Conv,ENV_trad.Fir_Conv.MC_trial),2),'s--','linewidth',2,'color',[i/size(ExpRes.bat_uratio_energy,1),0,0]);
            hold on;
            pic_fpp(i) = plot(ENV_fpp.Sweep.Conv_energy_sum/ENV_fpp.Sweep.Bat_energy_sum, ...
            100*mean(reshape(ExpRes_fpp.bat_uratio_energy(i,:,:),ENV_fpp.Sweep.Stat.Conv,ENV_fpp.Var_Conv.MC_trial),2),'o-.','linewidth',2,'color',[i/size(ExpRes.bat_uratio_energy,1),0,0]);
            xlabel('Normalized Aggregate Converter Rating');
            ylabel('Battery Capacity Utilization (%)');
            str1 = strcat({'Battery Capacity Variation'},{' '},...
                {num2str(100*ENV.Sweep.Bat{1}.qlim_var(i)/ENV.Sweep.Bat{1}.qlim_mu(i))},{'%'});
 %           str2 = {};
            title(str1,'Fontsize',9);
            legend([pic(i),pic_trad(i),pic_fpp(i)],'LSV-PPP','LSV-PPP2','LS-HiPPP','Fontsize',10);
            grid on;
            grid minor;
        end
    % xlabel('Sum of Converter energy ratio in Traditional DPP');
    % ylabel('Bat Utiliation');
 
        %% Plot Sytem Efficiency  
        figure(16);
        for i = 1:size(ExpRes.bat_uratio_energy,1)
            subplot(4,2,i);
            pic(i) = plot((ENV.Fir_Conv.Num*ENV.Sweep.Conv.e_lim)/(CL.Stat.Bat_num*CL.Bat{1}.curlim*CL.Bat{1}.volt),...
            100-100*mean(reshape(ExpRes.total_energy_process(i,:,:)*0.15./ExpRes.max_output_energy(i,:,:),ENV.Sweep.Stat.Conv,ENV.Fir_Conv.MC_trial),2),'d-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_energy,1),0,0]);
            hold on;
            pic_trad(i) = plot(ENV_trad.Sweep.Conv_energy_sum/ENV_trad.Sweep.Bat_energy_sum, ...
            100-100*mean(reshape(ExpRes_trad.total_energy_process(i,:,:)*0.15./ExpRes_trad.max_output_energy(i,:,:),ENV_trad.Sweep.Stat.Conv,ENV_trad.Fir_Conv.MC_trial),2),'s--','linewidth',2,'color',[i/size(ExpRes_trad.bat_uratio_energy,1),0,0]);
            hold on;
%           pic_fpp(i) = plot(ENV_trad.Sweep.Total_energy/Bat_energy_sum, ...
%           0.85*ones(1,8),'o-.','linewidth',2,'color',[i/size(ExpRes_fpp.bat_uratio_energy,1),0,0]);
            pic_fpp(i) = plot(ENV_fpp.Sweep.Conv_energy_sum/ENV_fpp.Sweep.Bat_energy_sum, ...
            100-100*mean(reshape(ExpRes_fpp.total_energy_process(i,:,:)*0.15./ExpRes_fpp.max_output_energy(i,:,:),ENV_fpp.Sweep.Stat.Conv,ENV_fpp.Var_Conv.MC_trial),2),'o-.','linewidth',2,'color',[i/size(ExpRes_fpp.bat_uratio_energy,1),0,0]);
            xlabel('Normalized Aggregate Converter Rating');
            ylabel('System Energy Efficiency (%)');
            str1 = strcat({'Battery Capacity Variation'},{' '},...
                {num2str(100*ENV.Sweep.Bat{1}.qlim_var(i)/ENV.Sweep.Bat{1}.qlim_mu(i))},{'%'});
            title(str1,'Fontsize',9);
            legend([pic(i),pic_trad(i),pic_fpp(i)],'LSV-PPP','LSV-PPP2','LS-HiPPP','Fontsize',10);
            grid on;
            grid minor;
        end
%     for i = 1:size(ExpRes.bat_uratio_energy,1)
%         pic(i) = plot((ENV.Var_Conv.Num*ENV.Sweep.Conv.p_lim + sum(ENV.Avg_Conv.p_lim_vec))/(CL.Stat.Bat_num*CL.Bat{1}.qlim*CL.Bat{1}.volt),...
%             1-mean(reshape(ExpRes.total_energy_process(i,:,:)*0.85./ExpRes.max_output_energy(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2),'d-','linewidth',2);
%         hold on;
%     end
 
    % xlabel('Sum of Converter energy Rate in Traditional DPP');
    % ylabel('Bat Utiliation');
%     xlabel('Converter Rating / Battery Rating');
%     ylabel('System energy Efficiency (Assume 85% Converter Efficiency)');
%     legend(pic,'5% Capacity Variance','10% Capacity Variance','15% Capacity Variance',...
%         '20% Capacity Variance','Fontsize',10);
%    legend(pic_trad,'Trad Bat energy Var = 0.01','Trad Bat energy Var = 0.02','Trad Bat energy Var = 0.03','Trad Bat energy Var = 0.04','Fontsize',10);
%    xlim([0 1.5]);
%     grid on;
%     grid minor;
    
%% Plot the utilzation at only 25% Variaion
figure();
        for i = 5:5
            pic_25 = plot((ENV.Fir_Conv.Num*ENV.Sweep.Conv.e_lim)/(CL.Stat.Bat_num*CL.Bat{1}.curlim*CL.Bat{1}.volt),...
            100*mean(reshape(ExpRes.bat_uratio_energy(i,:,:),ENV.Sweep.Stat.Conv,ENV.Fir_Conv.MC_trial),2),'d-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_energy,1),0,0]);
            hold on;
            pic_trad_25 = plot(ENV_trad.Sweep.Conv_energy_sum/ENV_trad.Sweep.Bat_energy_sum, ...
            100*mean(reshape(ExpRes_trad.bat_uratio_energy(i,:,:),ENV_trad.Sweep.Stat.Conv,ENV_trad.Fir_Conv.MC_trial),2),'s--','linewidth',2,'color',[i/size(ExpRes.bat_uratio_energy,1),0,0]);
            hold on;
            pic_fpp_25 = plot(ENV_fpp.Sweep.Conv_energy_sum/ENV_fpp.Sweep.Bat_energy_sum, ...
            100*mean(reshape(ExpRes_fpp.bat_uratio_energy(i,:,:),ENV_fpp.Sweep.Stat.Conv,ENV_fpp.Var_Conv.MC_trial),2),'o-.','linewidth',2,'color',[i/size(ExpRes.bat_uratio_energy,1),0,0]);
            xlabel('Normalized Aggregate Converter Rating');
            ylabel('Battery Energy Utilization (%)');
            title('Battery Capacity Variation 25%','Fontsize',9);
            legend([pic_25,pic_trad_25,pic_fpp_25],'LSV-PPP','LSV-PPP2','LS-HiPPP','Fontsize',10,'Location',[0.3,0.25,0.3,0.1]);
            grid on;
            grid minor;
        end
        ylim([70 100]);
 
%% Plot Sytem Efficiency at only 25% Variation
 figure();
        for i = 5:5
            pic_25 = plot((ENV.Fir_Conv.Num*ENV.Sweep.Conv.e_lim)/(CL.Stat.Bat_num*CL.Bat{1}.curlim*CL.Bat{1}.volt),...
            100-100*mean(reshape(ExpRes.total_energy_process(i,:,:)*0.15./ExpRes.max_output_energy(i,:,:),ENV.Sweep.Stat.Conv,ENV.Fir_Conv.MC_trial),2),'d-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_energy,1),0,0]);
            hold on;
            pic_trad_25 = plot(ENV_trad.Sweep.Conv_energy_sum/ENV_trad.Sweep.Bat_energy_sum, ...
            100-100*mean(reshape(ExpRes_trad.total_energy_process(i,:,:)*0.15./ExpRes_trad.max_output_energy(i,:,:),ENV_trad.Sweep.Stat.Conv,ENV_trad.Fir_Conv.MC_trial),2),'s--','linewidth',2,'color',[i/size(ExpRes_trad.bat_uratio_energy,1),0,0]);
            hold on;
%           pic_fpp(i) = plot(ENV_trad.Sweep.Total_energy/Bat_energy_sum, ...
%           0.85*ones(1,8),'o-.','linewidth',2,'color',[i/size(ExpRes_fpp.bat_uratio_energy,1),0,0]);
            pic_fpp_25 = plot(ENV_fpp.Sweep.Conv_energy_sum/ENV_fpp.Sweep.Bat_energy_sum, ...
            100-100*mean(reshape(ExpRes_fpp.total_energy_process(i,:,:)*0.15./ExpRes_fpp.max_output_energy(i,:,:),ENV_fpp.Sweep.Stat.Conv,ENV_fpp.Var_Conv.MC_trial),2),'o-.','linewidth',2,'color',[i/size(ExpRes_fpp.bat_uratio_energy,1),0,0]);
            xlabel('Normalized Aggregate Converter Rating');
            ylabel('System Energy Efficiency (%)');
            title('Battery Capacity Variation 25%','Fontsize',9);
            legend([pic_25,pic_trad_25,pic_fpp_25],'LSV-PPP','LSV-PPP2','LS-HiPPP','Fontsize',10);
            grid on;
            grid minor;
        end    
    
    %% Plot the Battery Utilization ratio and Converter Utilization ratio Tradeoff
%     figure(12);
%     for i = 1:size(ExpRes.bat_uratio_energy,1)
%         pic(i) = plot(mean(reshape(ExpRes.bat_uratio_energy_conv(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2),...
%             mean(reshape(ExpRes.bat_uratio_energy(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2),'d-','linewidth',2);
%         hold on;
%     end
%     xlabel('Converter Utiliation');
%     ylabel('Bat Utiliation');
%     legend(pic,'Bat energy Var = 0.01','Bat energy Var = 0.02','Bat energy Var = 0.03','Bat energy Var = 0.04','Fontsize',10);
%     grid on;
%     grid minor;
% 
%     %% Plot the Processed energy
%     figure(13);
%     for i = 1:size(ExpRes.process_energy,1)
%         pic(i) = plot(ENV.Sweep.Conv.p_lim, mean(reshape(ExpRes.process_energy(i,:,:), ENV.Sweep.Stat.Conv, ENV.Var_Conv.MC_trial),2),'s-','linewidth',2);
%         hold on;
%     end
%     xlabel('Var Converter energy');
%     ylabel('Processed energy');
%     legend(pic,'Bat energy Var = 0.01','Bat energy Var = 0.02','Bat energy Var = 0.03','Bat energy Var = 0.04','Fontsize',10);
%     grid on;
%     grid minor;
% 
%     %% Plot the Maximum Output energy
%     figure(14);
%     for i = 1:size(ExpRes.max_energy,1)
%         pic(i) = plot(ENV.Sweep.Conv.p_lim, mean(reshape(ExpRes.max_energy(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2),'*-','linewidth',2);
%         hold on;
%     end
%     xlabel('Var Converter energy');
%     ylabel('Max Output energy');
%     legend(pic,'Bat energy Var = 0.01','Bat energy Var = 0.02','Bat energy Var = 0.03','Bat energy Var = 0.04','Fontsize',10);
%     grid on;
%     grid minor;
% 
%     %% Plot the partial energy ratio
%     figure(15);
%     ppp_ratio = cell(size(ExpRes.max_energy,1),size(ExpRes.max_energy,2));
%     ppp_ratio_mean = [];
%     for k0 = 1:size(ExpRes.max_energy,1)
%         for k1 = 1:size(ExpRes.max_energy,2)
%             for k2 = 1:size(ExpRes.max_energy,3)
%                 if(ExpRes.max_energy(k0,k1,k2) ~= 0)
%                     ppp_ratio{k0,k1}(end+1) = ExpRes.process_energy(k0,k1,k2)./ExpRes.max_energy(k0,k1,k2);
%                 else
%                     ppp_ratio{k0,k1}(end+1) = 1;   % indicating we have to use full energy processing
%                 end
%                 ppp_ratio_mean(k0,k1) = mean(ppp_ratio{k0,k1});
%             end
%         end
%     end
%     for i = 1:size(ExpRes.max_energy,1)
%         pic(i) = plot(ENV.Sweep.Conv.p_lim, ppp_ratio_mean(i,:),'d--','linewidth',2);
%         hold on;
%     end
%     xlabel('Var Converter energy');
%     ylabel('Partial energy Ratio');
%     legend(pic,'Bat energy Var = 0.01','Bat energy Var = 0.02','Bat energy Var = 0.03','Bat energy Var = 0.04','Fontsize',10);
%     grid on;
%     grid minor;
end
