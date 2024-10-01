function y = SelfCheck_ENV6(CL,ENV6)
%UNTITLED Summary of this function goes here
for i = 1:CL.Stat.Bat_num
    if((length(ENV6.Sweep.Bat{i}.qlim_pop_var)~=ENV6.Sweep.Stat.BatPopVar)||...
        (length(ENV6.Sweep.Bat{i}.qlim_test_var)~=ENV6.Sweep.Stat.TstBatT))  
    fprintf('Assignment to ENV6.Sweep is wrong');
    pause;
    end
end
end

