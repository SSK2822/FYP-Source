function [Max_Drawdown0] = Calmar_Var_p_index(weeks, index_returns)
[~, Num_weeks] = size(index_returns);
wk_return = index_returns(:, Num_weeks/2+1:end);
Drawdown = zeros(1, weeks);

for week = 1:weeks
    rj = wk_return(:, 1:week);
    accu_return = zeros(1, week);
    for i = 1:week
        accu_return(i) = sum(rj(:, 1:i));
    end
    Drawdown(week) = max(accu_return)-sum(rj);
end
Max_Drawdown0 = max(Drawdown);
end
