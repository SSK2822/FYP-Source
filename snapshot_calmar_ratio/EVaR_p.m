function [EVaR0] = EVaR_p(weeks, wk_return_d1, xt_all, rho, theta)
%whole dataset input, x is decision matrix(multi-period)
%used for ex-post version
[~, N] = size(wk_return_d1);
wk_return = wk_return_d1(:, N/2+1:end);

exp_sum = 0;
for week = 1:weeks
    rj = wk_return(:, week); 
    exp_sum = exp_sum + exp(-rj'*xt_all(:, week)/rho);
end
EVaR0 = rho * log(1/(weeks)*exp_sum) - rho * log(1-theta) ;
end
