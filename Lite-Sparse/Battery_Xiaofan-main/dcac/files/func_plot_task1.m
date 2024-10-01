function y = func_plot_task1(CL,ENV,ExpRes,ENV_trad,ExpRes_trad)
%% Multiple Figures on Battery Utilization VS Aggregate Converter Rating, LS-HiPPP and C-PPP
a = figure(11);
set(a,'units','centimeters','position',[5,5,16,12]);
picture_num = 4;      % type 4 for 5%10%15%20%, type 8 for 5%10%15%20%25%50%75%100%
base_color = [0 0 1];
for i = 1:picture_num 
    subplot(picture_num/2,2,i);
    pic(i) = plot(ENV.Sweep.Conv_power_sum/ENV.Sweep.Bat_power_sum,...
        100*mean(reshape(ExpRes.bat_uratio_power(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2),'d-','linewidth',2,'color',i/picture_num*base_color);
    hold on;
    pic_trad(i) = plot(ENV_trad.Sweep.Conv_power_sum/ENV_trad.Sweep.Bat_power_sum, ...
    100*mean(reshape(ExpRes_trad.bat_uratio_power(i,:,:),ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),2),'s--','linewidth',2,'color',i/picture_num*base_color);
    xlabel('Normalized Aggregate Converter Rating');
    ylabel('Battery Power Utilization (%)');
    str1 = strcat({'Battery Power Capability'});
    str2 = strcat({'Variation'},{' '},...
        {num2str(100*ENV.Sweep.Bat{1}.curlim_var(i)/ENV.Sweep.Bat{1}.curlim_mu(i))},{'%'});
    str3 = {' '};
    title([str1; str2; str3],'Fontsize',9);
    legend([pic(i),pic_trad(i)],'LS-HiPPP','C-PPP','Fontsize',10,'Location','southeast');
    grid on;
    grid minor;
end

    % xlabel('Sum of Converter Power ratio in Traditional DPP');
    % ylabel('Bat Utiliation');
%% Multiple Figures on System Efficiency VS Aggregate Converter Rating, LS-HiPPP and C-PPP
a = figure();
set(a,'units','centimeters','position',[5,5,16,12]);
picture_num = 4;   % type 4 for 5%10%15%20%, type 8 for 5%10%15%20%25%50%75%100%
base_color = [0 1 0];
for i = 1:picture_num 
    subplot(picture_num/2,2,i);
    pic(i) = plot(ENV.Sweep.Conv_power_sum/ENV.Sweep.Bat_power_sum,...
        100-100*mean(reshape(ExpRes.total_power_process(i,:,:)*0.15./ExpRes.max_output_power(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2),'d-','linewidth',2,'color',i/picture_num*base_color);
    hold on;
    pic_trad(i) = plot(ENV_trad.Sweep.Conv_power_sum/ENV_trad.Sweep.Bat_power_sum, ...
    100-100*mean(reshape(ExpRes_trad.total_power_process(i,:,:)*0.15./ExpRes_trad.max_output_power(i,:,:),ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),2),'s--','linewidth',2,'color',i/picture_num*base_color);
    xlabel('Normalized Aggregate Converter Rating');
    ylabel('System Power Efficiency (%)');
%     str1 = strcat({'Battery Power Capability Variation'},{' '},...
%         {num2str(100*ENV.Sweep.Bat{1}.curlim_var(i)/ENV.Sweep.Bat{1}.curlim_mu(i))},{'%'});
    str1 = strcat({'Battery Power Efficiency'});
    str2 = strcat({'Variation'},{' '},...
        {num2str(100*ENV.Sweep.Bat{1}.curlim_var(i)/ENV.Sweep.Bat{1}.curlim_mu(i))},{'%'});
    str3 = {' '};
    title([str1; str2; str3],'Fontsize',9);
    legend([pic(i),pic_trad(i)],'LS-HiPPP','C-PPP','Fontsize',10,'Location','east');
    grid on;
    grid minor;
end
%% LS-HiPPP Figures on Converter Processed Power and Battery Utilization Tradeoff 
figure();
    for i = 1:4
        conv_uratio_vec = 100*mean(reshape(ExpRes.total_power_process(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2)./ENV.Sweep.Bat_power_sum;
        bat_uratio_vec = 100*mean(reshape(ExpRes.bat_uratio_power(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2);
        pic5(i) = plot(conv_uratio_vec(3:end),bat_uratio_vec(3:end),'d-','linewidth',2);
        hold on;
    end
    xlabel('Converter Processed Power (%)');
    ylabel('Bat Utilization (%)');
    legend(pic5,'5% Heterogeneity','10% Heterogeneity','15% Heterogeneity','20% Heterogeneity','Fontsize',10);
    title('LS-HiPPP');    
    grid on;
    grid minor;
    xlim([0 50]);
%% C-PPP Figures on Converter Processed Power and Battery Utilization Tradeoff 
figure();
    for i = 1:4
        pic5(i) = plot(100*mean(reshape(ExpRes_trad.total_power_process(i,:,:),ENV_trad.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2)./ENV.Sweep.Bat_power_sum...
            ,100*mean(reshape(ExpRes_trad.bat_uratio_power(i,:,:),ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),2),'d-','linewidth',2);
        hold on;
    end
    xlabel('Converter Processed Power (%)');
    ylabel('Bat Utilization (%)');
    legend(pic5,'5% Heterogeneity','10% Heterogeneity','15% Heterogeneity','20% Heterogeneity','Fontsize',10);
    title('C-PPP'); 
    grid on;
    grid minor;
    xlim([0 50]);
%% C-PPP Figures and LS-HiPPP uratio-processed power Tradeoff Comparison at 20% heterogeneity 
figure();
for i = 4:4
    conv_uratio_vec = 100*mean(reshape(ExpRes.total_power_process(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2)./ENV.Sweep.Bat_power_sum;
    bat_uratio_vec = 100*mean(reshape(ExpRes.bat_uratio_power(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2);
    pic_lshipppp = plot(conv_uratio_vec,bat_uratio_vec,'d-','linewidth',2);
    hold on;
end
for i = 4:4
    conv_uratio_trad_vec = 100*mean(reshape(ExpRes_trad.total_power_process(i,:,:),ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),2)./ENV_trad.Sweep.Bat_power_sum;
    bat_uratio_trad_vec = 100*mean(reshape(ExpRes_trad.bat_uratio_power(i,:,:),ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),2);
    pic_cpppp = plot(conv_uratio_trad_vec,bat_uratio_trad_vec ,'+--','linewidth',2);
    hold on;
end
    xlabel('Converter Processed Power (%)');
    ylabel('Battery Power Utilization (%)');
    legend([pic_lshipppp, pic_cpppp], 'LS-HiPPP', 'C-PPP');
    grid on;
    grid minor;
%     xlim([0 0.5]);

%% C-PPP Figures and LS-HiPPP system efficiecny - processed power Tradeoff Comparison at 20% heterogeneity 
figure();
for i = 4:4
    %conv_uratio_vec = 100*mean(reshape(ExpRes.total_power_process(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2)./ENV.Sweep.Bat_power_sum;
    sys_eff_vec =  100-100*mean(reshape(ExpRes.total_power_process(i,:,:)*0.15./ExpRes.max_output_power(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2);
    bat_uratio_vec = 100*mean(reshape(ExpRes.bat_uratio_power(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2);
    % pic_lshipppp = plot(sys_eff_vec,bat_uratio_vec,'s:','Color',[0.4940, 0.1840, 0.5560],'linewidth',2);
    pic_lshipppp = plot(bat_uratio_vec,sys_eff_vec,'s:','Color',[0.4940, 0.1840, 0.5560],'linewidth',2);
    hold on;
end
for i = 4:4
    %conv_uratio_trad_vec = 100*mean(reshape(ExpRes_trad.total_power_process(i,:,:),ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),2)./ENV_trad.Sweep.Bat_power_sum;
    sys_eff_trad_vec = 100-100*mean(reshape(ExpRes_trad.total_power_process(i,:,:)*0.15./ExpRes_trad.max_output_power(i,:,:),ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),2);
    bat_uratio_trad_vec = 100*mean(reshape(ExpRes_trad.bat_uratio_power(i,:,:),ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),2);
    % pic_cpppp = plot(sys_eff_trad_vec ,bat_uratio_trad_vec ,'o-.','Color',[0.4660, 0.6740, 0.1880],'linewidth',2);
    pic_cpppp = plot(bat_uratio_trad_vec,sys_eff_trad_vec,'o-.','Color',[0.4660, 0.6740, 0.1880],'linewidth',2);
    hold on;
end
    ylabel('System Power Efficiency (%)');
    xlabel('Battery Power Utilization (%)');
    legend([pic_lshipppp, pic_cpppp], 'LS-HiPPP', 'C-PPP');
    grid on;
    grid minor;
%     xlim([0 0.5])
end