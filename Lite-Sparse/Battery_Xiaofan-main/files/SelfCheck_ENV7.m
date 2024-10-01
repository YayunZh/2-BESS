function y = SelfCheck_ENV7(CL,ENV)
%UNTITLED Summary of this function goes here
for i = 1:CL.Stat.Bat_num
    if(length(ENV.Sweep.Bat{i}.qlim_pop_var)~=ENV.Sweep.Stat.BatPopVar)||...
        (length(ENV.Sweep.Bat{i}.qlim_diagnosis_var)~=ENV.Sweep.Stat.BatDiagnosis)
    printf('Assignment to ENV7.Sweep is wrong');
    pause;
    end
end

if(length(ENV.Sweep.DiscountFactor)~=ENV.Sweep.Stat.BatDiagnosis)
    printf('Assignment to ENV7.Sweep is wrong');
end

end
