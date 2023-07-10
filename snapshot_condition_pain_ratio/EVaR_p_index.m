function [EVaR0] = EVaR_p_index(weeks, wk_return_b1, rho, theta)
%whole dataset input, x is decision matrix(multi-period)
%used for ex-post version
[~, N] = size(wk_return_b1);
exp_sum = 0;
for week = N/2+1:weeks 
    rj = wk_return_b1(:, week); 
    exp_sum = exp_sum + exp(-rj/rho);
end
EVaR0 = rho * log(1/(weeks-1)*exp_sum) - rho * log(1-theta);
end
