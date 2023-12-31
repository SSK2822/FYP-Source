function [CDD_Drawdown0_a] = CDD_Var_a(weeks, wk_return_d1, x, rho, theta)
Drawdown = zeros(1, weeks-1); 

for week = 1:weeks-1
    rj = wk_return_d1(:, 1:week); 
    accu_return = zeros(1, week);
    for i = 1:week
        rkk = wk_return_d1(:, 1:i);
        accu_return(i) = sum(rkk'*x) ;
    end
    Drawdown(week) = max(accu_return)-sum(rj'*x);
end

maxi = max(Drawdown-rho,0);
val = (weeks-1)*(1-theta);
CDD_Drawdown0_a = rho + 1/(val)*sum(maxi);
end 