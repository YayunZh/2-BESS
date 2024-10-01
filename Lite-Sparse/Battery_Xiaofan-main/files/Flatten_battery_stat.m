function y = Flatten_battery_stat(CL)
    mu = CL.Bat{1}.qlim_mu;
    sigma = CL.Bat{1}.qlim_var;
    temp = norminv(linspace(1e-4,1-1e-4,CL.Stat.Bat_num+1),CL.Bat{1}.qlim_mu, CL.Bat{1}.qlim_var);   %Normal inverse cumulative distribution function
    for k = 1:CL.Stat.Bat_num
       x_int = linspace(temp(k),temp(k+1),100);
       y_int = normpdf(x_int,CL.Bat{1}.qlim_mu, CL.Bat{1}.qlim_var).*x_int;% p * x
       CL.Bat{k}.qlim = CL.Stat.Bat_num*trapz(x_int,y_int);% 9 * integral(p(x)x) - expected value
    end
    y = CL;
end