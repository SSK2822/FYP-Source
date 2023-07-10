function [CVaR0] = CVaR_p(weeks, wk_return_d1, xt_all, rho, theta)
%weeks: the ith weeks of test data 
%whole dataset input, x is decision matrix(multi-period)
%used for ex-post version
[~, N] = size(wk_return_d1);
wk_return = wk_return_d1(:, N/2+1:end);

max_sum = 0;
for week = 1:weeks
    rj = wk_return(:, week); 
    max_sum = max_sum + max((-1 * rj' * xt_all(:,week) - rho), 0);
end
% CVaR0 = rho * log(1/(weeks)*exp_sum) - rho * log(1-theta) ;
CVaR0 = rho + (1/(weeks*(1-theta)))*max_sum;
end