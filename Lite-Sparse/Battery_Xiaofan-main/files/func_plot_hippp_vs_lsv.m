function y = func_plot_hippp_vs_lsv(CL,ENV,ExpRes,ENV_trad,ExpRes_trad)
%% Plot the utilzation at only 25% Variaion
% sort: i-5,va-25, 6; i-6,va-50, start-; i-7,va-75;
start = [4,6,14,15];
% start = [1,1,1,1];
figure(1)
        for n = 4:7
            i = n-3;
            subplot(4,1,i);
            x = ENV.Sweep.Conv_energy_sum/(CL.Stat.Bat_num*CL.Bat{1}.qlim_mu*CL.Bat{1}.volt);
            x = x(start(i):end);
            y = 100*mean(reshape(ExpRes.bat_uratio_energy(i,:,:),ENV.Sweep.Stat.Conv,ENV.Fir_Conv.MC_trial),2);
            y = reshape(y(start(i):end),1,numel(y(start(i):end)));
            pic_lsv = plot(x,...
            y,'d-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_energy,1),0,0]);
            hold on;
            x_hi = ENV_trad.Sweep.Conv_energy_sum/ENV_trad.Sweep.Bat_energy_sum;
            x_hi = x_hi(start(i):end);
            y_hi = 100*mean(reshape(ExpRes_trad.bat_uratio_energy(i,:,:),ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),2);
            y_hi = reshape(y_hi(start(i):end),1,numel(y_hi(start(i):end)));
            pic_hippp = plot(x_hi, y_hi,'s--','linewidth',2,'color',[i/size(ExpRes.bat_uratio_energy,1),0,0]);
            xlabel('Normalized Aggregate Converter Rating');
            ylabel('Battery Energy Utilization (%)');
            str1 = strcat({'Battery Capacity Variation'},{' '},...
                {num2str(100*ENV.Sweep.Bat{1}.qlim_var(n)/ENV.Sweep.Bat{1}.qlim_mu(n))},{'%'});
 %           str2 = {};
            title(str1,'Fontsize',9);            
            legend([pic_lsv,pic_hippp],'LSV-PPP','LS-HiPPP','Fontsize',10);
            grid on;
            grid minor;
            xlim([0,1]);
        end
        % ylim([70 100]);
 
%% Plot Sytem Efficiency at only 25% Variation
 figure(2);
        for n = 4:7
            i = n - 3;
            subplot(4,1,i)
            x = ENV.Sweep.Conv_energy_sum/(CL.Stat.Bat_num*CL.Bat{1}.qlim_mu*CL.Bat{1}.volt);
            x = x(start(i):end);
            y = 100-100*mean(reshape(ExpRes.total_energy_process(i,:,:)*0.15./ExpRes.max_output_energy(i,:,:),ENV.Sweep.Stat.Conv,ENV.Fir_Conv.MC_trial),2);
            y = reshape(y(start(i):end),1,numel(y(start(i):end)));
            pic_lsv = plot(x,...
            y,'d-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_energy,1),0,0]);
            hold on;
            x_hi = ENV_trad.Sweep.Conv_energy_sum/ENV_trad.Sweep.Bat_energy_sum;
            x_hi = x_hi(start(i):end);
            y_hi = 100-100*mean(reshape(ExpRes_trad.total_energy_process(i,:,:)*0.15./ExpRes_trad.max_output_energy(i,:,:),ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),2);
            y_hi = reshape(y_hi(start(i):end),1,numel(y_hi(start(i):end)));
            pic_hippp = plot(x_hi, ...
            y_hi,'s--','linewidth',2,'color',[i/size(ExpRes_trad.bat_uratio_energy,1),0,0]);
            xlabel('Normalized Aggregate Converter Rating');
            ylabel('System Energy Efficiency (%)');
            str1 = strcat({'Battery Capacity Variation'},{' '},...
                {num2str(100*ENV.Sweep.Bat{1}.qlim_var(n)/ENV.Sweep.Bat{1}.qlim_mu(n))},{'%'});
            title(str1,'Fontsize',9);            
            legend([pic_lsv,pic_hippp],'LSV-PPP','LS-HiPPP','Fontsize',10);
            xlim([0,1]);
            grid on;
            grid minor;
        end    
 