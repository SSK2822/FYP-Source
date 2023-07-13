function [martin_drawdown_index] = Martin_Var_p_index(week, index_returns)
% Calculate Martin ratio for index data
[~, num_weeks] = size(index_returns);
index_returns_test = index_returns(:, num_weeks/2+1:end);
drawdown = zeros(1, week); 

for j = 1:week
    rj = index_returns_test(:, 1:j);
    cumulative_return = zeros(1, j);
    for i = 1:j
        cumulative_return(i) = sum(rj(:, 1:i));
    end
    drawdown(j) = (max(cumulative_return)-sum(rj))^2;
end
martin_drawdown_index = sum(drawdown)/week;
end
