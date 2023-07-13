function [average_drawdown] = Pain_Var_p(week, asset_returns, portfolio_weights_all)
set_returns);
asset_returns_test = asset_returns(:, num_weeks/2+1:end);
drawdown = zeros(1, week); 

for j = 1:week
    rj = asset_returns_test(:, 1:j);
    cumulative_return = zeros(1, j);
    for i = 1:j
        rkk = asset_returns_test(:, 1:i);
        xtt = portfolio_weights_all(:, 1:i);
        cumulative_return(i) =  sum(diag(rkk'*xtt));
    end
    
    drawdown(j) = max(cumulative_return)-sum(diag(rj'*portfolio_weights_all(:, 1:week))); % week period's maximum drawndown
end
average_drawdown = sum(drawdown)/week;

end
