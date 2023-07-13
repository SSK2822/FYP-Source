function [Max_Drawdown0] = Calmar_Var_p(weeks, assets_return, portfolio_all)
[~, N] = size(assets_return);
wk_return = assets_return(:, N/2+1:end);
Drawdown = zeros(1, weeks); 

for week = 1:weeks
    rj = wk_return(:, 1:week);
    accu_return = zeros(1, week);
    for i = 1:week
        rkk = wk_return(:, 1:i);
        xtt = portfolio_all(:, 1:i);
        accu_return(i) =  sum(diag(rkk'*xtt));
    end
    
    Drawdown(week) = max(accu_return)-sum(diag(rj'*portfolio_all(:, 1:week))); %week期的maximum drawndown
end
Max_Drawdown0 = max(Drawdown);

end

