function y = SelfCheck_ENV1_2(CL,ENV1_2)
%UNTITLED Summary of this function goes here
for i = 1:CL.Stat.Bat_num
    if((length(ENV1_2.Sweep.Bat{i}.curlim_var)~=ENV1_2.Sweep.Stat.Bat))  
        fprintf('Assignment to ENV1_2.Sweep is wrong');
        pause;
    end
end
if((length(ENV1_2.Sweep.Conv.p_lim)~=ENV1_2.Sweep.Stat.Conv))  
    fprintf('Assignment to ENV1_2.Sweep is wrong');
    pause;
end
end

