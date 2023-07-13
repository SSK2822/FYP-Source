function [CDD_Drawdown0] = CDD_Var_p_index(weeks, wk_return_b1, rho, theta)
[~, N] = size(wk_return_b1);
wk_return = wk_return_b1(:, N/2+1:end);
Drawdown = zeros(1, weeks);

for week = 1:weeks
    rj = wk_return(:, 1:week);
    accu_return = zeros(1, week);
    for i = 1:week
        accu_return(i) = sum(rj(:,1:i));
    end
    Drawdown(week) = max(accu_return)-sum(rj);
end
CDD_Drawdown0 = rho + 1/(weeks*(1-theta))*sum(max(Drawdown-rho,0));
end
