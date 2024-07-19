function Opt_PQ = func_PQ(SYS,TEST)

%% 
for a = 1 : TEST.Sweep.Stat.Bat % for different heterogeneous degree
    SYS.Bat{1}.curlim_var = TEST.Sweep.Bat{1}.curlim_var(a);
    SYS.Bat{1}.curlim_mu = TEST.Sweep.Bat{1}.curlim_mu(a);
    SYS = Flatten_battery_stat_cur(SYS); % calculate expected values

    for b = 1 : TEST.Sweep.Stat.Conv % for different normalized aggregate converter rating
        TEST.Conv.p_lim_singlevar = TEST.Sweep.Conv.p_lim(b);

        for c = 1:TEST.Sweep.Stat.phi
            TEST.t_delay = TEST.Sweep.phi(c);

            OptRes(a,b,c) = func_energyflow_PQ(SYS, TEST);
 
            Power_capability(a,b,c) =  OptRes(a,b,c).I_max;

        end
    end
end
        
Opt_PQ.power_capability = Power_capability;
Opt_PQ.OptRes = OptRes;

save('PQI_max.mat','Power_capability');
save('PQoutput.mat','OptRes');

end

%%
function y = Flatten_battery_stat_cur(SYS)
    temp = norminv(linspace(1e-4,1-1e-4,SYS.Stat.Bat_num+1),SYS.Bat{1}.curlim_mu, SYS.Bat{1}.curlim_var);   %Normal inverse cumulative distribution function
    for k = 1:SYS.Stat.Bat_num
       x_int = linspace(temp(k),temp(k+1),100);
       y_int = normpdf(x_int,SYS.Bat{1}.curlim_mu, SYS.Bat{1}.curlim_var).*x_int;% p * x
       SYS.Bat{k}.curlim = SYS.Stat.Bat_num*trapz(x_int,y_int);% 9 * integral(p(x)x) - expected value
       SYS.Bat{k}.qlim = SYS.Bat{k}.curlim*10;
    end
    y = SYS;
end

