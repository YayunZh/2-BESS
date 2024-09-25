function Opt_MC = func_MC(SYS,TEST)

load('QL0.1.mat');
load('Conv_Fir_mat0.1.mat')

for k0 = 1:4 %: TEST.Sweep.Stat.Bat       % for different heterogeneous degree
    [SYS, Bat_Info_MC] = Update_k0_loop(k0,SYS,TEST);

    for k1 = 1 : TEST.Sweep.Stat.Conv   % for different normalized aggregate converter rating
        TEST.Conv.p_lim_singlevar = TEST.Sweep.Conv.p_lim(k1);
        QL_connection = QL{k0,k1};
        P_Fir_mat = Conv_Fir_mat{k0,k1};

        parfor k2 = 1:TEST.Conv.MC_trial % for MC simu

            OptResMC(k0,k1,k2) = func_energyflow_MC(Bat_Info_MC(:,k2),SYS,TEST,QL_connection,P_Fir_mat);

            Up(k0,k1,k2) = OptResMC(k0,k1,k2).up;
            Eff(k0,k1,k2) = OptResMC(k0,k1,k2).efficience;
            Eff2(k0,k1,k2) = OptResMC(k0,k1,k2).efficience2;
        end 

        Up_avg(k0,k1) = sum(Up(k0,k1,:))/TEST.Conv.MC_trial;
        Eff_avg(k0,k1) = sum(Eff(k0,k1,:))/TEST.Conv.MC_trial;
        Eff2_avg(k0,k1) = sum(Eff2(k0,k1,:))/TEST.Conv.MC_trial;
    end
end

Opt_MC.power_Utilization = Up_avg;
Opt_MC.efficience = Eff_avg;
Opt_MC.efficience2 = Eff2_avg;
Opt_MC.MC = OptResMC;

end

%% MC simu information generate
function [SYS,Bat_Info_MC] = Update_k0_loop(k0,SYS,TEST) 
    for i = 1:SYS.Stat.Bat_num
        SYS.Bat{i}.curlim_var = TEST.Sweep.Bat{i}.curlim_var(k0);
        SYS.Bat{i}.curlim_mu = TEST.Sweep.Bat{i}.curlim_mu(k0);

        pd = makedist('Normal','mu',TEST.Sweep.Bat{i}.curlim_mu(k0),'sigma',TEST.Sweep.Bat{i}.curlim_var(k0));
        t = truncate(pd,0,inf);
        Bat_Info_MC(i,:) = random(t,1,TEST.Conv.MC_trial);
    end
    for j = 1:TEST.Conv.MC_trial % sort batteries curlim
        Bat_Info_MC(:,j) = sort(Bat_Info_MC(:,j));
    end
end


