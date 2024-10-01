function y = func_plot_task1_4_power(CL,ENV,ExpRes,ENV_trad,ExpRes_trad,ENV_fpp,ExpRes_fpp)
    figure();
    picture_num = 4;   % type 4 for 5%10%15%20%, type 8 for 5%10%15%20%25%50%75%100%
    base_color = [0 0 1];   % base color is blue
        i = 4;
    pic(i) = plot(ENV.Sweep.Conv_power_sum/ENV.Sweep.Bat_power_sum,...
        100*mean(reshape(ExpRes.bat_uratio_power(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2),'d-','linewidth',2,'color',i/picture_num*base_color);
    hold on;
    pic_trad(i) = plot(ENV_trad.Sweep.Conv_power_sum/ENV_trad.Sweep.Bat_power_sum, ...
    100*mean(reshape(ExpRes_trad.bat_uratio_power(i,:,:),ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),2),'s--','linewidth',2,'color',i/picture_num*base_color);
    hold on;
    pic_fpp(i) = plot(ENV_fpp.Sweep.Conv_power_sum/ENV_fpp.Sweep.Bat_power_sum, ...
    100*mean(reshape(ExpRes_fpp.bat_uratio_power(i,:,:),ENV_fpp.Sweep.Stat.Conv,ENV_fpp.Var_Conv.MC_trial),2),'o-.','linewidth',2,'color',i/picture_num*base_color);
    xlabel('Normalized Aggregate Converter Rating');
    ylabel('Battery Power Utilization (%)');
%     str1 = strcat({'Battery Power Capability Variation'},{' '},...
%         {num2str(100*ENV.Sweep.Bat{1}.curlim_var(i)/ENV.Sweep.Bat{1}.curlim_mu(i))},{'%'});
%           str2 = {};
%    title(str1,'Fontsize',9);
    legend([pic(i),pic_trad(i),pic_fpp(i)],'LS-HiPPP','C-PPP','Full Processing','Fontsize',10);
    grid on;
    grid minor;
    ylim([70 100]);
end

