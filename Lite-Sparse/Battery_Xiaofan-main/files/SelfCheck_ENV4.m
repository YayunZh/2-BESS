function y = SelfCheck_ENV4(CL,ENV4)
%UNTITLED Summary of this function goes here
for i = 1:CL.Stat.Bat_num
    if((length(ENV4.Sweep.Bat{i}.qlim_var)~=ENV4.Sweep.Stat.Bat))  
        fprintf('Assignment to ENV4.Sweep is wrong');
        pause;
    end
end
if((length(ENV4.Sweep.Conv.e_lim)~=ENV4.Sweep.Stat.Conv))  
    fprintf('Assignment to ENV4.Sweep is wrong');
    pause;
end
end

