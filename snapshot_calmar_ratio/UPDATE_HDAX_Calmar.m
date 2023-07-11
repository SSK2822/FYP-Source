% Clear workspace and command window
clc,clear

% Set display format to long
format long;

% Load data from file
asset_prices = load('wk_price_HDAX', ',')'; % Asset prices matrix (939 days of data)
index_prices = load('b1_HDAX.txt', ','); % Index prices
output_file = './UPDATE_HDAX_Calmar/'; % Output file path

risk_free_rate = load('rf_00_17.txt',','); % Risk-free rate data

% Convert prices to returns
[asset_returns, ~] = price2ret(asset_prices', [], 'Periodic'); % Simple return
asset_returns = asset_returns';
[index_returns, ~] = price2ret(index_prices, [], 'Periodic'); % Simple return

% Get number of assets and weeks
[num_assets, num_weeks] = size(asset_returns);

% Variable initial
theta = 0.95;
portfolio_weights = zeros(num_assets, 1);
portfolio_weights_all = zeros(num_assets, num_weeks/2);
equal_weights = 1 / num_assets * ones(num_assets, 1);
calmar_ratio_a_all = zeros(1, num_weeks/2);
calmar_ratio_p_all = zeros(1, num_weeks/2);
nolinear_all = zeros(1, num_weeks/2);
rho_all = zeros(1, num_weeks/2);
ksi_all = zeros(num_assets, num_weeks/2);
portfolio_return_all = zeros(1,num_weeks/2);

% Loop over the second half of weeks
for week = (num_weeks/2+1):num_weeks 
    % Display current week
    disp(week)
    
    % Get training and test data for asset returns
    asset_returns_train = asset_returns(:, 1:week-1); 
    asset_returns_test = asset_returns(:, week); 
    current_risk_free_rate = risk_free_rate(week); % Risk-free rate
    mean_return_train = mean(asset_returns_train, 2); % Mean return
    
    if (mean_return_train < current_risk_free_rate) 
        fprintf('do not trade in week %d\n', week)
        portfolio_weights = zeros(num_assets, 1); 
    else 
        [portfolio_weights, ksi] = NNN_IG_linear(week, asset_returns, risk_free_rate, mean_return_train); 
    end     
    portfolio_weights_all(:, week-num_weeks/2) = portfolio_weights;
    
    % Calmar ratio for test data
    portfolio_return_test = portfolio_weights' * asset_returns_test - current_risk_free_rate;
    portfolio_return_all(week-num_weeks/2) = portfolio_return_test;
    
    max_drawdown_test = Calmar_Var_p(week-num_weeks/2, asset_returns, portfolio_weights_all);
    calmar_ratio_test = sum(portfolio_return_all) / (week-num_weeks/2) / max_drawdown_test;
    calmar_ratio_p_all(week-num_weeks/2) = calmar_ratio_test;
    
    % Calmar ratio for training data
    portfolio_return_train = portfolio_weights' * mean(asset_returns_train, 2) - current_risk_free_rate;
    calmar_ratio_train = portfolio_return_train / Calmar_Var_a(week, asset_returns, portfolio_weights);
    calmar_ratio_a_all(week-num_weeks/2) = calmar_ratio_train;
    

    % Save the results of portfolio weights to a file
    fid0 = fopen([output_file, 'xt_', num2str(week), '.txt'], 'w'); % Output file name
    for ii = 1:num_assets
        fprintf(fid0, '%.10f\t', portfolio_weights(ii, 1));
        fprintf(fid0, '\r\n');
    end
    fclose(fid0);
    
end
% Calmar ratio for updated dataset
calmar_ratio_update = calmar_ratio_p_all* (num_weeks/18)^0.5; 
