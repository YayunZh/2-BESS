function y = SelfCheck_ENV9(CL,ENV)
%UNTITLED Summary of this function goes here
for i = 1:CL.Stat.Bat_num
    if(length(ENV.Sweep.Var_conv.lambda_scaling)~=ENV.Sweep.Stat.LambdaScaling)||...
        (length(ENV.Sweep.Bat{i}.qlim_diagnosis_var)~=ENV.Sweep.Stat.BatDiagnosis)
    fprintf('Assignment to ENV9.Sweep is wrong');
    pause;
    end
end

if(length(ENV.Sweep.DiscountFactor)~=ENV.Sweep.Stat.BatDiagnosis)
    fprintf('Assignment to ENV9.Sweep is wrong');
    pause;
end

end
