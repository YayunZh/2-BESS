function Opt_connection = func_connection_design(SYS,TEST)

%% 
for a = 1 : TEST.Sweep.Stat.Bat % for different heterogeneous degree
    SYS.Bat{1}.curlim_var = TEST.Sweep.Bat{1}.curlim_var(a);
    SYS.Bat{1}.curlim_mu = TEST.Sweep.Bat{1}.curlim_mu(a);
    SYS = Flatten_battery_stat_cur(SYS); % calculate expected values

    for b = 1 : TEST.Sweep.Stat.Conv % for different normalized aggregate converter rating
        TEST.Conv.p_lim_singlevar = TEST.Sweep.Conv.p_lim(b);

        OptRes(a,b) = func_energyflow_connection(SYS, TEST);
        QL{a,b} = OptRes(a,b).Bat_use;
        Conv_Fir_mat{a,b} = OptRes(a,b).Fir_mat ;

    end
end

save('QL.mat','QL');
save('Conv_Fir_mat.mat','Conv_Fir_mat');

Opt_connection.OptRes = OptRes;
Opt_connection.QL = QL;
Opt_connection.Conv_fir_mat = Conv_Fir_mat;

end


%%
function y = Flatten_battery_stat_cur(SYS)
    temp = norminv(linspace(1e-4,1-1e-4,SYS.Stat.Bat_num+1),SYS.Bat{1}.curlim_mu, SYS.Bat{1}.curlim_var);   % Normal inverse cumulative distribution function
    for k = 1:SYS.Stat.Bat_num
       x_int = linspace(temp(k),temp(k+1),100);
       y_int = normpdf(x_int,SYS.Bat{1}.curlim_mu, SYS.Bat{1}.curlim_var).*x_int;
       SYS.Bat{k}.curlim = SYS.Stat.Bat_num*trapz(x_int,y_int);
       SYS.Bat{k}.capalim = SYS.Bat{k}.curlim*10;
    end
    y = SYS;
end



