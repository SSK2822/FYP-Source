function [c,ceq] = nonlinear_constrain(u, week, asset_returns, theta, risk_free_rate)
% asset_returns: the whole dataset, risk_free_rate: vector
% DSR used for ex-ante version
[num_assets, ~] = size(asset_returns);
y = u(1:num_assets,:);
ksi = u(num_assets+2:end, :);
ita = u(num_assets+1,:);

c = zeros(1, 2*(week-1));
con_sum = zeros(1, week-1);

for j = 1:week-1
    rj = asset_returns(:, 1:j); 
    con_sum(j) = ksi(j)-sum(rj'*y);
    c(j) = con_sum(j)-1;
    c(week-1+j) = 0-con_sum(j);
end
% c = [con_sum-1; 0-con_sum];
ceq = [];

end
