function y = Update_Bat_stat(Bat_Info,CL)
    for l = 1:CL.Stat.Bat_num
        %CL.Bat{l}.qlim = mvnrnd([CL.Bat{l}.qlim_mu,CL.Bat{l}.volt_mu], diag([CL.Bat{l}.qlim_var,CL.Bat{l}.volt_var]),1);
         CL.Bat{l}.qlim = Bat_Info(l);
    end
    CL.Cap{3}.volt = CL.Bus{1}.volt - CL.Bat{1}.volt - CL.Bat{2}.volt - CL.Bat{3}.volt...
    - CL.Bat{4}.volt - CL.Bat{5}.volt - CL.Bat{6}.volt...
    - CL.Bat{7}.volt - CL.Bat{8}.volt - CL.Bat{9}.volt;     % Decide the capactior votlages
    y = CL;
end