function y = func_plot_task1_4(CL,ENV,ExpRes,ENV_trad,ExpRes_trad,ENV_fpp,ExpRes_fpp)
    %% Multiple Plots on the Utilization ratio VS aggregate converter rating -- LS-HiPPP, C-PPP and F-PPP 
    a = figure();
    set(a,'units','centimeters','position',[5,5,16,12]);
    picture_num = 4;   % type 4 for 5%10%15%20%, type 8 for 5%10%15%20%25%50%75%100%
    base_color = [0 0 1];   % base color is blue
        for i = 1:picture_num
            subplot(picture_num/2,2,i);
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
            ylim([70 100]);
            str1 = strcat({'Battery Power Capability'});
            str2 = strcat({'Variation'},{' '},...
                {num2str(100*ENV.Sweep.Bat{1}.curlim_var(i)/ENV.Sweep.Bat{1}.curlim_mu(i))},{'%'});
            str3 = {' '};
            title([str1,str2,str3],'Fontsize',9);
            legend([pic(i),pic_trad(i),pic_fpp(i)],'LS-HiPPP','C-PPP','Full Processing','Fontsize',10,'Location','southwest');
            grid on;
            grid minor;
        end

    % xlabel('Sum of Converter Power ratio in Traditional DPP');
    % ylabel('Bat Utiliation');
 
      %% Multiple Plots on the system efficiency VS aggregate converter rating -- LS-HiPPP, C-PPP and F-PPP   
        figure();
        picture_num = 4; % type 4 for 5%10%15%20%, type 8 for 5%10%15%20%25%50%75%100%
        base_color = [0 1 0];
        for i = 1:picture_num
            subplot(picture_num/2,2,i);
            pic(i) = plot(ENV.Sweep.Conv_power_sum/ENV.Sweep.Bat_power_sum,...
                100-100*mean(reshape(ExpRes.total_power_process(i,:,:)*0.15./(ExpRes.max_output_power(i,:,:))...
                ,ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2),'d-','linewidth',2,'color',i/picture_num*base_color);
            ylim([70 100]);
            hold on;
            pic_trad(i) = plot(ENV_trad.Sweep.Conv_power_sum/ENV_trad.Sweep.Bat_power_sum, ...
                100-100*mean(reshape(ExpRes_trad.total_power_process(i,:,:)*0.15./(ExpRes_trad.max_output_power(i,:,:))...
                ,ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),2),'s--','linewidth',2,'color',i/picture_num*base_color);
            ylim([70 100]);
            hold on;
%           pic_fpp(i) = plot(ENV_trad.Sweep.Total_power/Bat_power_sum, ...
%           0.85*ones(1,8),'o-.','linewidth',2,'color',[i/size(ExpRes_fpp.bat_uratio_power,1),0,0]);
            pic_fpp(i) = plot(ENV_fpp.Sweep.Conv_power_sum/ENV_fpp.Sweep.Bat_power_sum, ...
                100-100*mean(reshape(ExpRes_fpp.total_power_process(i,:,:)*0.15./(ExpRes_fpp.max_output_power(i,:,:))...
                ,ENV_fpp.Sweep.Stat.Conv,ENV_fpp.Var_Conv.MC_trial),2),'o-.','linewidth',2,'color',i/picture_num*base_color);
            ylim([70 100]);
            xlabel('Normalized Aggregate Converter Rating');
            ylabel('System Power Efficiency (%)');
            str1 = strcat({'Battery Power Capability Variation'},{' '},...
                {num2str(100*ENV.Sweep.Bat{1}.curlim_var(i)/ENV.Sweep.Bat{1}.curlim_mu(i))},{'%'});
            title(str1,'Fontsize',9);
            legend([pic(i),pic_trad(i),pic_fpp(i)],'LS-HiPPP','C-PPP','Full Processing','Fontsize',10);
            grid on;
            grid minor;
            
        end    

        
      %% Plot for the system efficiency VS aggregate converter rating at 20% uncertianty
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
            str1 = strcat({'Battery Power Capability Variation'},{' '},...
                {num2str(100*ENV.Sweep.Bat{1}.curlim_var(i)/ENV.Sweep.Bat{1}.curlim_mu(i))},{'%'});
 %           str2 = {};
            title(str1,'Fontsize',9);
            legend([pic(i),pic_trad(i),pic_fpp(i)],'LS-HiPPP','C-PPP','Full Processing','Fontsize',10);
            grid on;
            grid minor;
    % xlabel('Sum of Converter Power ratio in Traditional DPP');
    % ylabel('Bat Utiliation');

        
 %% Plot Processed Power VS Aggregate Converter Rating
 figure();
for i = 4:4
    pic_25 = plot(ENV.Sweep.Conv_power_sum/ENV.Sweep.Bat_power_sum,...
    100*mean(reshape(ExpRes.total_power_process(i,:,:),ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2)/ENV.Sweep.Bat_power_sum,'d-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_power,1),0,0]);
    hold on;
    pic_trad_25 = plot(ENV_trad.Sweep.Conv_power_sum/ENV_trad.Sweep.Bat_power_sum, ...
    100*mean(reshape(ExpRes_trad.total_power_process(i,:,:),ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),2)/ENV.Sweep.Bat_power_sum,'s--','linewidth',2,'color',[i/size(ExpRes_trad.bat_uratio_power,1),0,0]);
    hold on;
%           pic_fpp(i) = plot(ENV_trad.Sweep.Total_power/Bat_power_sum, ...
%           0.85*ones(1,8),'o-.','linewidth',2,'color',[i/size(ExpRes_fpp.bat_uratio_power,1),0,0]);
    pic_fpp_25 = plot(ENV_fpp.Sweep.Conv_power_sum/ENV_fpp.Sweep.Bat_power_sum, ...
    100*mean(reshape(ExpRes_fpp.total_power_process(i,:,:),ENV_fpp.Sweep.Stat.Conv,ENV_fpp.Var_Conv.MC_trial),2)/ENV.Sweep.Bat_power_sum,'o-.','linewidth',2,'color',[i/size(ExpRes_fpp.bat_uratio_power,1),0,0]);
    xlabel('Normalized Aggregate Converter Rating');
    ylabel('Normalized Processed Power (%)');
    title('Battery Power Capability Variation 20%','Fontsize',9);
    legend([pic_25,pic_trad_25,pic_fpp_25],'LS-HiPPP','C-PPP','Full Processing','Fontsize',10);
    grid on;
    grid minor;
end  
        
%% Plot the Utilzation VS battery heterogeneity at only 20% Converter Power Rating (8th data)
% 1 assume r (conv raitng/battery raintg)= 0.2, utilization   VS heterogenity 
figure();
bat_rating_ind = 8; 
ENV.Sweep.Conv_power_sum(bat_rating_ind)/ENV.Sweep.Bat_power_sum
        for i = 1:4  % use ENV.Sweep.Stat.Bat if want to plot the data for all heterogeneity
            
            uratio2batheter = 100*mean(reshape(ExpRes.bat_uratio_power(i,:,:),... 
            ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2);
            uratio2batheter_max = 100*prctile(reshape(ExpRes.bat_uratio_power(i,:,:),... 
            ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),90,2);
            uratio2batheter_min = 100*prctile(reshape(ExpRes.bat_uratio_power(i,:,:),... 
            ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),10,2);
            heter(i) = 100*ENV.Sweep.Bat{1}.curlim_var(i)/ENV.Sweep.Bat{1}.curlim_mu(i);
            u2heter(i) = uratio2batheter(bat_rating_ind);
            u2heter_max(i) = uratio2batheter_max(bat_rating_ind);
            u2heter_min(i) = uratio2batheter_min(bat_rating_ind);
            %pic_25 = plot(heter, u2heter, 'd-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_power,1),0,0]);
            hold on;
            
            uratio2batheter_trad = 100*mean(reshape(ExpRes_trad.bat_uratio_power(i,:,:),... 
            ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),2);
            uratio2batheter_trad_max = 100*prctile(reshape(ExpRes_trad.bat_uratio_power(i,:,:),... 
            ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),90,2);
            uratio2batheter_trad_min = 100*prctile(reshape(ExpRes_trad.bat_uratio_power(i,:,:),... 
            ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),10,2);
            heter(i) = 100*ENV_trad.Sweep.Bat{1}.curlim_var(i)/ENV_trad.Sweep.Bat{1}.curlim_mu(i);
            u2heter_trad(i) = uratio2batheter_trad(bat_rating_ind);
            u2heter_trad_max(i) = uratio2batheter_trad_max(bat_rating_ind);
            u2heter_trad_min(i) = uratio2batheter_trad_min(bat_rating_ind);  
            %pic_trad_25 = plot(heter, u2heter_trad, 'd-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_power,1),0,0]);

            uratio2batheter_fpp = 100*mean(reshape(ExpRes_fpp.bat_uratio_power(i,:,:),... 
            ENV_fpp.Sweep.Stat.Conv,ENV_fpp.Var_Conv.MC_trial),2);
            uratio2batheter_fpp_max = 100*prctile(reshape(ExpRes_fpp.bat_uratio_power(i,:,:),... 
            ENV_fpp.Sweep.Stat.Conv,ENV_fpp.Var_Conv.MC_trial),90,2);
            uratio2batheter_fpp_min = 100*prctile(reshape(ExpRes_fpp.bat_uratio_power(i,:,:),... 
            ENV_fpp.Sweep.Stat.Conv,ENV_fpp.Var_Conv.MC_trial),10,2);
            heter(i) = 100*ENV_fpp.Sweep.Bat{1}.curlim_var(i)/ENV_fpp.Sweep.Bat{1}.curlim_mu(i);
            u2heter_fpp(i) = uratio2batheter_fpp(bat_rating_ind);
            u2heter_fpp_max(i) = uratio2batheter_fpp_max(bat_rating_ind);
            u2heter_fpp_min(i) = uratio2batheter_fpp_min(bat_rating_ind);
            %pic_fpp_25 = plot(heter, u2heter_fpp, 'd-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_power,1),0,0]);
        end
        pic_25 =  errorbar(heter, u2heter, u2heter_min - u2heter, u2heter_max - u2heter,'o-','linewidth',2,'color',[0.9290, 0.6940, 0.1250]);
%         pic_25 = plot(heter, u2heter, );
        hold on;
        pic_trad_25 = errorbar(heter, u2heter_trad, u2heter_trad_min - u2heter_trad, u2heter_trad_max - u2heter_trad, 'x:','linewidth',2,'color',[0.4940, 0.1840, 0.5560]);
%       pic_trad_25 = plot(heter, u2heter_trad, 'x:','linewidth',2,'color','r');
        % hold on;
        % pic_fpp_25 = errorbar(heter, u2heter_fpp, u2heter_fpp_min - u2heter_fpp, u2heter_fpp_max - u2heter_fpp, 'x-.','linewidth',2,'color',[0.4660, 0.6740, 0.1880]);
%       pic_fpp_25 = plot(heter, u2heter_fpp, '^-.','linewidth',2,'color','g');
        
        xlabel('Battery Power Heterogeneity (%)');
        ylabel('Battery Power Utilization (%)');
        title('Converter Rating = 20% of Battery Rating','Fontsize',9);
        % legend([pic_25,pic_trad_25,pic_fpp_25],'LS-HiPPP','C-PPP','Full Processing','Fontsize',10,'Location',[0.3,0.25,0.3,0.1]);
         legend([pic_25,pic_trad_25],'LS-HiPPP','C-PPP','Fontsize',10,'Location',[0.3,0.25,0.3,0.1]);
        grid on;
        grid minor;
        xlim([5 20]);

%% Plot the Pout VS battery heterogeneity at only 20% Converter Power Rating (8th data)
% 1 assume r (conv raitng/battery raintg)= 0.2, utilization   VS heterogenity 
figure();
bat_rating_ind = 8;
        for i = 1:size(ExpRes.bat_uratio_power,1)
            
            uratio2batheter = 100*mean(reshape(ExpRes.max_output_power(i,:,:),... 
            ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2)/ENV.Sweep.Bat_power_sum;
            uratio2batheter_max = 100*prctile(reshape(ExpRes.max_output_power(i,:,:),... 
            ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),90,2)/ENV.Sweep.Bat_power_sum;
            uratio2batheter_min = 100*prctile(reshape(ExpRes.max_output_power(i,:,:),... 
            ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),10,2)/ENV.Sweep.Bat_power_sum;
            heter(i) = 100*ENV.Sweep.Bat{1}.curlim_var(i)/ENV.Sweep.Bat{1}.curlim_mu(i);
            u2heter(i) = uratio2batheter(bat_rating_ind);
            u2heter_max(i) = uratio2batheter_max(bat_rating_ind);
            u2heter_min(i) = uratio2batheter_min(bat_rating_ind);
            %pic_25 = plot(heter, u2heter, 'd-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_power,1),0,0]);          

            uratio2batheter_trad = 100*mean(reshape(ExpRes_trad.max_output_power(i,:,:),... 
            ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),2)/ENV_trad.Sweep.Bat_power_sum;
            uratio2batheter_trad_max = 100*prctile(reshape(ExpRes_trad.max_output_power(i,:,:),... 
            ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),90,2)/ENV_trad.Sweep.Bat_power_sum;
            uratio2batheter_trad_min = 100*prctile(reshape(ExpRes_trad.max_output_power(i,:,:),... 
            ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),10,2)/ENV_trad.Sweep.Bat_power_sum;
            heter(i) = 100*ENV_trad.Sweep.Bat{1}.curlim_var(i)/ENV_trad.Sweep.Bat{1}.curlim_mu(i);
            u2heter_trad(i) = uratio2batheter_trad(bat_rating_ind);
            u2heter_trad_max(i) = uratio2batheter_trad_max(bat_rating_ind);
            u2heter_trad_min(i) = uratio2batheter_trad_min(bat_rating_ind);  
            %pic_trad_25 = plot(heter, u2heter_trad, 'd-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_power,1),0,0]);

            uratio2batheter_fpp = 100*mean(reshape(ExpRes_fpp.max_output_power(i,:,:),... 
            ENV_fpp.Sweep.Stat.Conv,ENV_fpp.Var_Conv.MC_trial),2)/ENV_fpp.Sweep.Bat_power_sum;
            uratio2batheter_fpp_max = 100*prctile(reshape(ExpRes_fpp.max_output_power(i,:,:),... 
            ENV_fpp.Sweep.Stat.Conv,ENV_fpp.Var_Conv.MC_trial),90,2)/ENV_fpp.Sweep.Bat_power_sum;
            uratio2batheter_fpp_min = 100*prctile(reshape(ExpRes_fpp.max_output_power(i,:,:),... 
            ENV_fpp.Sweep.Stat.Conv,ENV_fpp.Var_Conv.MC_trial),10,2)/ENV_fpp.Sweep.Bat_power_sum;
            heter(i) = 100*ENV_fpp.Sweep.Bat{1}.curlim_var(i)/ENV_fpp.Sweep.Bat{1}.curlim_mu(i);
            u2heter_fpp(i) = uratio2batheter_fpp(bat_rating_ind);
            u2heter_fpp_max(i) = uratio2batheter_fpp_max(bat_rating_ind);
            u2heter_fpp_min(i) = uratio2batheter_fpp_min(bat_rating_ind);
            %pic_fpp_25 = plot(heter, u2heter_fpp, 'd-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_power,1),0,0]);
        end
        pic_25 =  errorbar(heter, u2heter, u2heter_min - u2heter, u2heter_max - u2heter,'o-','linewidth',2,'color',[0.9290, 0.6940, 0.1250]);
%         pic_25 = plot(heter, u2heter, );
        hold on;
        pic_trad_25 = errorbar(heter, u2heter_trad, u2heter_trad_min - u2heter_trad, u2heter_trad_max - u2heter_trad, 'x:','linewidth',2,'color',[0.4940, 0.1840, 0.5560]);
%       pic_trad_25 = plot(heter, u2heter_trad, 'x:','linewidth',2,'color','r');
%        hold on;
%        pic_fpp_25 = errorbar(heter, u2heter_fpp, u2heter_fpp_min - u2heter_fpp, u2heter_fpp_max - u2heter_fpp, 'x-.','linewidth',2,'color',[0.4660, 0.6740, 0.1880]);
%       pic_fpp_25 = plot(heter, u2heter_fpp, '^-.','linewidth',2,'color','g');
        
        xlabel('Battery Power Heterogeneity (%)');
        ylabel('Maximum Power Output (%)');
        title('Converter Rating = 20% of Battery Rating','Fontsize',9);
        % legend([pic_25,pic_trad_25,pic_fpp_25],'LS-HiPPP','C-PPP','Full Processing','Fontsize',10,'Location',[0.3,0.25,0.3,0.1]);
         legend([pic_25,pic_trad_25],'LS-HiPPP','C-PPP','Fontsize',10,'Location',[0.3,0.25,0.3,0.1]);
        grid on;
        grid minor;
        xlim([5 20]);

%% Plot the Pout VS battery heterogeneity at only 100% Converter Power Rating (20th data)
% 1 assume r (conv raitng/battery raintg)= 0.2, utilization   VS heterogenity 
figure();
bat_rating_ind = 4;
ENV.Sweep.Bat{1}.curlim_var(bat_rating_ind)
        for i = 1:size(ExpRes.bat_uratio_power,1)
            
            uratio2batheter = 100*mean(reshape(ExpRes.max_output_power(i,:,:),... 
            ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2)/ENV.Sweep.Bat_power_sum;
            uratio2batheter_max = 100*prctile(reshape(ExpRes.max_output_power(i,:,:),... 
            ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),90,2)/ENV.Sweep.Bat_power_sum;
            uratio2batheter_min = 100*prctile(reshape(ExpRes.max_output_power(i,:,:),... 
            ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),10,2)/ENV.Sweep.Bat_power_sum;
            heter(i) = 100*ENV.Sweep.Bat{1}.curlim_var(i)/ENV.Sweep.Bat{1}.curlim_mu(i);
            u2heter(i) = uratio2batheter(bat_rating_ind);
            u2heter_max(i) = uratio2batheter_max(bat_rating_ind);
            u2heter_min(i) = uratio2batheter_min(bat_rating_ind);
            %pic_25 = plot(heter, u2heter, 'd-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_power,1),0,0]);          

            uratio2batheter_trad = 100*mean(reshape(ExpRes_trad.max_output_power(i,:,:),... 
            ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),2)/ENV_trad.Sweep.Bat_power_sum;
            uratio2batheter_trad_max = 100*prctile(reshape(ExpRes_trad.max_output_power(i,:,:),... 
            ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),90,2)/ENV_trad.Sweep.Bat_power_sum;
            uratio2batheter_trad_min = 100*prctile(reshape(ExpRes_trad.max_output_power(i,:,:),... 
            ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),10,2)/ENV_trad.Sweep.Bat_power_sum;
            heter(i) = 100*ENV_trad.Sweep.Bat{1}.curlim_var(i)/ENV_trad.Sweep.Bat{1}.curlim_mu(i);
            u2heter_trad(i) = uratio2batheter_trad(bat_rating_ind);
            u2heter_trad_max(i) = uratio2batheter_trad_max(bat_rating_ind);
            u2heter_trad_min(i) = uratio2batheter_trad_min(bat_rating_ind);  
            %pic_trad_25 = plot(heter, u2heter_trad, 'd-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_power,1),0,0]);

            uratio2batheter_fpp = 100*mean(reshape(ExpRes_fpp.max_output_power(i,:,:),... 
            ENV_fpp.Sweep.Stat.Conv,ENV_fpp.Var_Conv.MC_trial),2)/ENV_fpp.Sweep.Bat_power_sum;
            uratio2batheter_fpp_max = 100*prctile(reshape(ExpRes_fpp.max_output_power(i,:,:),... 
            ENV_fpp.Sweep.Stat.Conv,ENV_fpp.Var_Conv.MC_trial),90,2)/ENV_fpp.Sweep.Bat_power_sum;
            uratio2batheter_fpp_min = 100*prctile(reshape(ExpRes_fpp.max_output_power(i,:,:),... 
            ENV_fpp.Sweep.Stat.Conv,ENV_fpp.Var_Conv.MC_trial),10,2)/ENV_fpp.Sweep.Bat_power_sum;
            heter(i) = 100*ENV_fpp.Sweep.Bat{1}.curlim_var(i)/ENV_fpp.Sweep.Bat{1}.curlim_mu(i);
            u2heter_fpp(i) = uratio2batheter_fpp(bat_rating_ind);
            u2heter_fpp_max(i) = uratio2batheter_fpp_max(bat_rating_ind);
            u2heter_fpp_min(i) = uratio2batheter_fpp_min(bat_rating_ind);
            %pic_fpp_25 = plot(heter, u2heter_fpp, 'd-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_power,1),0,0]);
        end
        pic_25 =  errorbar(heter, u2heter, u2heter_min - u2heter, u2heter_max - u2heter,'o-','linewidth',2,'color',[0.9290, 0.6940, 0.1250]);
%         pic_25 = plot(heter, u2heter, );
        hold on;
        pic_trad_25 = errorbar(heter, u2heter_trad, u2heter_trad_min - u2heter_trad, u2heter_trad_max - u2heter_trad, 'x:','linewidth',2,'color',[0.4940, 0.1840, 0.5560]);
%       pic_trad_25 = plot(heter, u2heter_trad, 'x:','linewidth',2,'color','r');
        hold on;
        pic_fpp_25 = errorbar(heter, u2heter_fpp, u2heter_fpp_min - u2heter_fpp, u2heter_fpp_max - u2heter_fpp, 'x-.','linewidth',2,'color',[0.4660, 0.6740, 0.1880]);
%       pic_fpp_25 = plot(heter, u2heter_fpp, '^-.','linewidth',2,'color','g');
        
        xlabel('Battery Power Heterogeneity (%)');
        ylabel('Maximum Power Output (%)');
        title('Converter Rating = 100% of Battery Rating','Fontsize',9);
        legend([pic_25,pic_trad_25,pic_fpp_25],'LS-HiPPP','C-PPP','Full Processing','Fontsize',10,'Location',[0.3,0.25,0.3,0.1]);
        % legend([pic_25,pic_trad_25],'LS-HiPPP','C-PPP','Fontsize',10,'Location',[0.3,0.25,0.3,0.1]);
        grid on;
        grid minor;
        xlim([5 20]);

        
%% Plot the System Efficiency VS battery heterogeneity at only 20% Converter Power Rating (8th data)
figure();
        for i = 1:size(ExpRes.bat_uratio_power,1)
            
            uratio2batheter = 100-100*mean(reshape(ExpRes.total_power_process(i,:,:)*0.15./ExpRes.max_output_power(i,:,:),...
            ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2);
            uratio2batheter_max = 100-100*prctile(reshape(ExpRes.total_power_process(i,:,:)*0.15./ExpRes.max_output_power(i,:,:),...
            ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),10,2);
            uratio2batheter_min = 100-100*prctile(reshape(ExpRes.total_power_process(i,:,:)*0.15./ExpRes.max_output_power(i,:,:),...
            ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),90,2);
            heter(i) = 100*ENV.Sweep.Bat{1}.curlim_var(i)/ENV.Sweep.Bat{1}.curlim_mu(i);
            u2heter(i) = uratio2batheter(8);
            u2heter_max(i) = uratio2batheter_max(8);
            u2heter_min(i) = uratio2batheter_min(8);
            %pic_25 = plot(heter, u2heter, 'd-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_power,1),0,0]);
            
            uratio2batheter_trad = 100-100*mean(reshape(ExpRes_trad.total_power_process(i,:,:)*0.15./ExpRes_trad.max_output_power(i,:,:),...
            ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),2);
            uratio2batheter_trad_max = 100-100*prctile(reshape(ExpRes_trad.total_power_process(i,:,:)*0.15./ExpRes_trad.max_output_power(i,:,:),...
            ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),10,2);
            uratio2batheter_trad_min = 100-100*prctile(reshape(ExpRes_trad.total_power_process(i,:,:)*0.15./ExpRes_trad.max_output_power(i,:,:),...
            ENV_trad.Sweep.Stat.Conv,ENV_trad.Var_Conv.MC_trial),90,2);
            heter(i) = 100*ENV_trad.Sweep.Bat{1}.curlim_var(i)/ENV_trad.Sweep.Bat{1}.curlim_mu(i);
            u2heter_trad(i) = uratio2batheter_trad(8);
            u2heter_trad_max(i) = uratio2batheter_trad_max(8);
            u2heter_trad_min(i) = uratio2batheter_trad_min(8);  
            %pic_trad_25 = plot(heter, u2heter_trad, 'd-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_power,1),0,0]);

            uratio2batheter_fpp = 100-100*mean(reshape(ExpRes_fpp.total_power_process(i,:,:)*0.15./ExpRes_fpp.max_output_power(i,:,:),...
            ENV_fpp.Sweep.Stat.Conv,ENV_fpp.Var_Conv.MC_trial),2);
            uratio2batheter_fpp_max = 100-100*prctile(reshape(ExpRes_fpp.total_power_process(i,:,:)*0.15./ExpRes_fpp.max_output_power(i,:,:),...
            ENV_fpp.Sweep.Stat.Conv,ENV_fpp.Var_Conv.MC_trial),10,2);
            uratio2batheter_fpp_min = 100-100*prctile(reshape(ExpRes_fpp.total_power_process(i,:,:)*0.15./ExpRes_fpp.max_output_power(i,:,:),...
            ENV_fpp.Sweep.Stat.Conv,ENV_fpp.Var_Conv.MC_trial),90,2);
            heter(i) = 100*ENV_fpp.Sweep.Bat{1}.curlim_var(i)/ENV_fpp.Sweep.Bat{1}.curlim_mu(i);
            u2heter_fpp(i) = uratio2batheter_fpp(8);
            u2heter_fpp_max(i) = uratio2batheter_fpp_max(8);
            u2heter_fpp_min(i) = uratio2batheter_fpp_min(8);
            %pic_fpp_25 = plot(heter, u2heter_fpp, 'd-','linewidth',2,'color',[i/size(ExpRes.bat_uratio_power,1),0,0]);
        end
        pic_25 =  errorbar(heter, u2heter, u2heter_min - u2heter, u2heter_max - u2heter,'o-','linewidth',2,'color',[0.9290, 0.6940, 0.1250]);
%         pic_25 = plot(heter, u2heter, );
        hold on;
        pic_trad_25 = errorbar(heter, u2heter_trad, u2heter_trad_min - u2heter_trad, u2heter_trad_max - u2heter_trad, 'x:','linewidth',2,'color',[0.4940, 0.1840, 0.5560]);
%       pic_trad_25 = plot(heter, u2heter_trad, 'x:','linewidth',2,'color','r');
        % hold on;
%         pic_fpp_25 = errorbar(heter, u2heter_fpp, u2heter_fpp_min - u2heter_fpp, u2heter_fpp_max - u2heter_fpp, 'x-.','linewidth',2,'color',[0.4660, 0.6740, 0.1880]);
%       pic_fpp_25 = plot(heter, u2heter_fpp, '^-.','linewidth',2,'color','g');
        
        xlabel('Battery Power Heterogeneity (%)');
        ylabel('System Power Efficiency (%)');
        title('Converter Rating = 20% of Battery Rating','Fontsize',9);
        % legend([pic_25,pic_trad_25,pic_fpp_25],'LS-HiPPP','C-PPP','Full Processing','Fontsize',10,'Location',[0.3,0.25,0.3,0.1]);
         legend([pic_25,pic_trad_25],'LS-HiPPP','C-PPP','Fontsize',10,'Location',[0.3,0.25,0.3,0.1]);
        grid on;
        grid minor;
        xlim([5 20]);
        
%%  utilization VS lambda under heterogeneity 5% 10% 15% 20%
figure();

 for i = 1:4
    lambda = (ENV.Sweep.Conv_power_sum - ENV.Sweep.Conv_power_sum_primary)/ENV.Sweep.Conv_power_sum_primary;
    uratio2lambda = 100*mean(reshape(ExpRes.bat_uratio_power(i,:,:),... 
    ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2);

    pic_25(i) = plot(lambda, uratio2lambda, '^-','linewidth',2,'color',[0,i/4,0]);
    hold on;  
 end
xlabel('Hierachy Factor (\lambda_H)');
ylabel('Battery Power Utilization (%)');
%        title('Utilization VS Hierachical','Fontsize',9);
legend(pic_25,'5% Heterogeneity','10% Heterogeneity','15% Heterogeneity', '20% Heterogeneity','Fontsize',10,'Location',[0.3,0.25,0.3,0.1]);
grid on;
grid minor;
xlim([0, 2.5])

%%  system efficiency VS lambda under heterogeneity 5% 10% 15% 20%
figure();

for i = 1:4 
    lambda = (ENV.Sweep.Conv_power_sum - ENV.Sweep.Conv_power_sum_primary)/ENV.Sweep.Conv_power_sum_primary;
    powereff2lambda = 100-100*mean(reshape(ExpRes.total_power_process(i,:,:)*0.15./ExpRes.max_output_power(i,:,:),...
        ENV.Sweep.Stat.Conv,ENV.Var_Conv.MC_trial),2);
    pic_25(i) = plot(lambda, powereff2lambda, 'v--','linewidth',2,'color',[0,0,i/4]);
    hold on;  
end
xlabel('Hierachy Factor (\lambda_H)');
ylabel('System Power Efficiency (%)');
%        title('Utilization VS Hierachical','Fontsize',9);
legend(pic_25,'5% Heterogeneity','10% Heterogeneity','15% Heterogeneity', '20% Heterogeneity','Fontsize',10,'Location',[0.3,0.25,0.3,0.1]);
grid on;
grid minor;
xlim([0, 2.5])

end

