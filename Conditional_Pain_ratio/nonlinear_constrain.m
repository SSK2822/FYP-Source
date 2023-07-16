function [c,ceq] = nonlinear_constrain(u, week, asset_returns, theta, risk_free_rate)
% asset_returns: the whole dataset, risk_free_rate: vector
% DSR used for ex-ante version
[num_assets, ~] = size(asset_returns);
y = u(1:num_assets,:);
ksi = u(num_assets+3:end, :);
rho = u(num_assets+2,:);
ita = u(num_assets+1,:);

con_sum = zeros(1, week-1);

for j = 1:week-1
    rj = asset_returns(:, 1:j); 
    con_sum(j) = max(ksi(j)-sum(rj'*y)-rho,0);
end
c = rho + 1/((week-1)*(1-theta))*sum(con_sum) - 1;
ceq = [];

end
