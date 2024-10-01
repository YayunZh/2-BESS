function y = func_plot_task1_4_short(CL,ENV,ExpRes,ENV_trad,ExpRes_trad,ENV_fpp,ExpRes_fpp)        
%% Plot the Utilzation VS battery heterogeneity at only 20% Converter Power Rating (8th data)
% 1 assume r (conv raitng/battery raintg)= 0.2, utilization   VS heterogenity 
figure();
bat_rating_ind = 4;
ENV.Sweep.Bat{1}.curlim_var(bat_rating_ind)
        for i = 1:ENV.Sweep.Stat.Conv
            
            uratio2batheter = 100*mean(reshape(ExpRes.bat_uratio_power(:,i,:),... 
            ENV.Sweep.Stat.Bat,ENV.Var_Conv.MC_trial),2);
            uratio2batheter_max = 100*prctile(reshape(ExpRes.bat_uratio_power(:,i,:),... 
            ENV.Sweep.Stat.Bat,ENV.Var_Conv.MC_trial),90,2);
            uratio2batheter_min = 100*prctile(reshape(ExpRes.bat_uratio_power(:,i,:),... 
            ENV.Sweep.Stat.Bat,ENV.Var_Conv.MC_trial),10,2);
            conv_rate(i) = 100*ENV.Sweep.Conv_power_sum(i)/ENV.Sweep.Bat_power_sum;
            u2conv(i) = uratio2batheter(bat_rating_ind);
            u2conv_max(i) = uratio2batheter_max(bat_rating_ind);
            u2conv_min(i) = uratio2batheter_min(bat_rating_ind);
            %pic_25 = plot(heter, u2conv, 'd-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_power,1),0,0]);
            hold on;
            
            uratio2batheter_trad = 100*mean(reshape(ExpRes_trad.bat_uratio_power(:,i,:),... 
            ENV_trad.Sweep.Stat.Bat,ENV_trad.Var_Conv.MC_trial),2);
            uratio2batheter_trad_max = 100*prctile(reshape(ExpRes_trad.bat_uratio_power(:,i,:),... 
            ENV_trad.Sweep.Stat.Bat,ENV_trad.Var_Conv.MC_trial),90,2);
            uratio2batheter_trad_min = 100*prctile(reshape(ExpRes_trad.bat_uratio_power(:,i,:),... 
            ENV_trad.Sweep.Stat.Bat,ENV_trad.Var_Conv.MC_trial),10,2);
            conv_rate(i) = 100*ENV.Sweep.Conv_power_sum(i)/ENV.Sweep.Bat_power_sum;
            u2conv_trad(i) = uratio2batheter_trad(bat_rating_ind);
            u2conv_trad_max(i) = uratio2batheter_trad_max(bat_rating_ind);
            u2conv_trad_min(i) = uratio2batheter_trad_min(bat_rating_ind);  
            %pic_trad_25 = plot(heter, u2conv_trad, 'd-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_power,1),0,0]);

            uratio2batheter_fpp = 100*mean(reshape(ExpRes_fpp.bat_uratio_power(:,i,:),... 
            ENV_fpp.Sweep.Stat.Bat,ENV_fpp.Var_Conv.MC_trial),2);
            uratio2batheter_fpp_max = 100*prctile(reshape(ExpRes_fpp.bat_uratio_power(:,i,:),... 
            ENV_fpp.Sweep.Stat.Bat,ENV_fpp.Var_Conv.MC_trial),90,2);
            uratio2batheter_fpp_min = 100*prctile(reshape(ExpRes_fpp.bat_uratio_power(:,i,:),... 
            ENV_fpp.Sweep.Stat.Bat,ENV_fpp.Var_Conv.MC_trial),10,2);
            conv_rate(i) = 100*ENV.Sweep.Conv_power_sum(i)/ENV.Sweep.Bat_power_sum;
            u2conv_fpp(i) = uratio2batheter_fpp(bat_rating_ind);
            u2conv_fpp_max(i) = uratio2batheter_fpp_max(bat_rating_ind);
            u2conv_fpp_min(i) = uratio2batheter_fpp_min(bat_rating_ind);
            %pic_fpp_25 = plot(heter, u2conv_fpp, 'd-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_power,1),0,0]);
        end
        pic_25 =  errorbar(conv_rate, u2conv, u2conv_min - u2conv, u2conv_max - u2conv,'o-','linewidth',2,'color',[0.9290, 0.6940, 0.1250]);
%         pic_25 = plot(heter, u2conv, );
        hold on;
        pic_trad_25 = errorbar(conv_rate, u2conv_trad, u2conv_trad_min - u2conv_trad, u2conv_trad_max - u2conv_trad, 'x:','linewidth',2,'color',[0.4940, 0.1840, 0.5560]);
%       pic_trad_25 = plot(heter, u2conv_trad, 'x:','linewidth',2,'color','r');
        % hold on;
        pic_fpp_25 = errorbar(conv_rate, u2conv_fpp, u2conv_fpp_min - u2conv_fpp, u2conv_fpp_max - u2conv_fpp, 'x-.','linewidth',2,'color',[0.4660, 0.6740, 0.1880]);
        % pic_fpp_25 = plot(heter, u2conv_fpp, '^-.','linewidth',2,'color','g');
        
        xlabel('Normalized Aggregate Converter Rating (%)');
        ylabel('Battery Power Utilization (%)');
        title('Battery Heterogeneity = 20%','Fontsize',9);
        legend([pic_25,pic_trad_25,pic_fpp_25],'LS-HiPPP','C-PPP','Full Processing','Fontsize',10,'Location',[0.3,0.25,0.3,0.1]);
        % legend([pic_25,pic_trad_25],'LS-HiPPP','C-PPP','Fontsize',10,'Location',[0.3,0.25,0.3,0.1]);
        grid on;
        grid minor;
%        xlim([5 20]);
%         ylim([70 100]);



end

