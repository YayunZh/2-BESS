function y = SelfCheck_ENV2(CL,ENV2)
%UNTITLED Summary of this function goes here
for i = 1:CL.Stat.Bat_num
    if((length(ENV2.Sweep.Bat{i}.curlim_pop_var)~=ENV2.Sweep.Stat.BatPopVar)||...
        (length(ENV2.Sweep.Bat{i}.curlim_test_var)~=ENV2.Sweep.Stat.TstBatT))  
    printf('Assignment to ENV2.Sweep is wrong');
    pause;
    end
end
end

