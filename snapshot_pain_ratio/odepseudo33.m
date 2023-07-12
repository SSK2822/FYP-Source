function objective = odepseudo33(u, mean_return, risk_free_rate, week, asset_returns_train)
% fmincon

% Get number of assets
[num_assets, ~] = size(asset_returns_train);

% Extract variables from input vector
y = u(1:num_assets , :);
ita = u(num_assets+1 , :);

% Define objective function
objective = -(mean_return'*y-risk_free_rate(week)*ita); 

end
