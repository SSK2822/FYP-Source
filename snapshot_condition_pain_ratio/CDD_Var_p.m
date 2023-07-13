function [CDD_Drawdown0] = CDD_Var_p(weeks, assets_return, portfolio_all, rho, theta)
[~, num_weeks] = size(assets_return);
week_return = assets_return(:, num_weeks/2+1:end);
Drawdown = zeros(1, weeks); 

for week = 1:weeks
    rj = week_return(:, 1:week);
    accu_return = zeros(1, week);
    for i = 1:week
        returnkk = week_return(:, 1:i);
        portfolio = portfolio_all(:, 1:i);
        accu_return(i) =  sum(diag(returnkk'*portfolio));
    end
    
    Drawdown(week) = max(accu_return)-sum(diag(rj'*portfolio_all(:, 1:week))); 
end

maxi = max(Drawdown-rho,0);
val = weeks*(1-theta);
CDD_Drawdown0 = rho + 1/(val)*sum(maxi);
end

