function [average_drawdown] = Pain_Var_p_index(week, index_returns)
[~, num_weeks] = size(index_returns);
index_returns_test = index_returns(:, num_weeks/2+1:end);
drawdown = zeros(1, week);

for j = 1:week
    rj = index_returns_test(:, 1:j); 
    cumulative_return = zeros(1, j);
    for i = 1:j
        cumulative_return(i) = sum(rj(:, 1:i));
    end
    drawdown(j) = max(cumulative_return)-sum(rj);
end
average_drawdown = sum(drawdown)/week;
end
