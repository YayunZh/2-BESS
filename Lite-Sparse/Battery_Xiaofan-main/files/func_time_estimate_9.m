function y = func_time_estimate_9(ENV,core)
    count_avg_design = ENV.Sweep.Stat.LambdaScaling*ENV.LambdaScaling.MC_trial*ENV.Avg_Conv.trial_num;
    count_var_design = ENV.Sweep.Stat.LambdaScaling*ENV.LambdaScaling.MC_trial*ENV.Sweep.Stat.BatDiagnosis*ENV.Diagnosis_Bat.MC_trial;
    unit_lp_time = 1;  % single LP takes 1 second 
    y = unit_lp_time*(count_avg_design + count_var_design)/core;
end

