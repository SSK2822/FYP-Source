function [CDD_Drawdown0] = CDD_Var_p(weeks, wk_return_d1, xt_all, rho, theta)
%whole dataset input, x is decision matrix(multi-period)
%used for ex-post version
%weeks: the ith weeks of test data
%计算后验CDD
[~, N] = size(wk_return_d1);
wk_return = wk_return_d1(:, N/2+1:end);
Drawdown = zeros(1, weeks); 

for week = 1:weeks
    rj = wk_return(:, 1:week);
    accu_return = zeros(1, week);
    for i = 1:week
        rkk = wk_return(:, 1:i);
        xtt = xt_all(:, 1:i);
        accu_return(i) =  sum(diag(rkk'*xtt));
    end
    
    Drawdown(week) = max(accu_return)-sum(diag(rj'*xt_all(:, 1:week))); %week期的maximum drawndown
end
CDD_Drawdown0 = rho + 1/(weeks*(1-theta))*sum(max(Drawdown-rho,0));
end
