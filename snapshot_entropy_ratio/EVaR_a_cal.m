function [EVaR0] = EVaR_a_cal(weeks, wk_return_d1, x, theta)
%wk_return_d1:the whole dataset, rb:vector
[~, EVaR0, a] = fmincon(@(rho) EVaR_a(weeks, wk_return_d1, x, rho, theta), ...
    0.01, [], [], [], [], 0, [], [], optimoptions('fmincon', 'Algorithm', 'sqp')...
            );
if(a ~= 1)
    disp(['a_cal:' num2str(a)]);
end
end 