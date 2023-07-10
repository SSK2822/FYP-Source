function [EVaR0] = EVaR_a(weeks, wk_return_d1, x, rho, theta)
%wk_return_d1:the whole dataset, rb:vector
%DSR used for ex-ante version
exp_sum = 0;

for week = 1:weeks-1
    rj = wk_return_d1(:, week); 
    exp_sum = exp_sum + exp(-rj'*x/rho);
end
EVaR0 = rho * log(1/(weeks-1)*exp_sum) - rho * log(1-theta);

end 