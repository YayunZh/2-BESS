function y = InputCheck(CL,ENV)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% if (size(ENV.Sweep.Bat.mu,2)~=size(ENV.Sweep.Bat.var,2))
%     fprintf('Error: The Battery Mu Dimension and Var Dimension mismatch \n');
%     pause;
% end
if(length(ENV.Avg_Conv.p_lim_vec) ~= ENV.Avg_Conv.Num)
    fprintf('Error: the Avg Converter number does not match the Number of Avg Converter Power Limit');
    pause;
end
end

