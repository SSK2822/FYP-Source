% Clear workspace and command window
clc,clear

% Set display format to long
format long;

% Get weekly prices of assets and index
asset_prices = load('wk_price_HDAX', ',')';
index_prices = load('b1_HDAX.txt', ',');

% Set output file path
output_file = './UPDATE_HDAX_Pain/';
output_file_rho = './UPDATE_HDAX_Pain_rho/';

% Load risk-free rate data
risk_free_rate = load('rf_00_17.txt',',');

% Convert prices to returns
[asset_returns, ~] = price2ret(asset_prices', [], 'Periodic');
asset_returns = asset_returns';
[index_returns, ~] = price2ret(index_prices, [], 'Periodic');

% Get number of assets and weeks
[num_assets, num_weeks] = size(asset_returns);

% Set parameter theta
theta = 0.95;

% Initialize variables
portfolio_weights = zeros(num_assets, 1);
portfolio_weights_all = zeros(num_assets, num_weeks/2);
equal_weights = 1 / num_assets * ones(num_assets, 1);
CDD_ratio_a_all = zeros(1, num_weeks/2);
CDD_ratio_p_all = zeros(1, num_weeks/2);
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
    
    % Get risk-free rate for current week
    current_risk_free_rate = risk_free_rate(week);
    
    % Get mean return of assets for training data
    mean_return_train = mean(asset_returns_train, 2);
    
    % Open output file for writing portfolio weights
    fid0 = fopen([output_file, 'xt_', num2str(week), '.txt'], 'w');
    
    % Check if any asset has mean return higher than risk-free rate
    if (mean_return_train < current_risk_free_rate)
        % If not, do not trade and set portfolio weights to zero
        fprintf('do not trade in week %d\n', week)
        portfolio_weights = zeros(num_assets, 1); 
    else 
        % If yes, use NNN_IG_linear function to get portfolio weights and ksi
        [portfolio_weights, rho, nolinear] = NNN_IG(week, asset_returns, risk_free_rate, mean_return_train);
     end     
    
    % Save portfolio weights for current week
    portfolio_weights_all(:, week-num_weeks/2) = portfolio_weights;
    rho_all(week-N/2) = rho;

    % Write portfolio weights to output file
    for ii = 1:num_assets
        fprintf(fid0, '%.10f\t', portfolio_weights(ii, 1));
        fprintf(fid0, '\r\n');
    end
    
    % Close output file
    fclose(fid0); 
    
    fid0 = fopen([output_file_rho, 'rho_', num2str(week), '.txt'], 'w');
    fprintf(fid0, '%.10f\t', rho);
    fclose(fid0);

    % Calculate portfolio return for test data
    portfolio_return_test = portfolio_weights' * asset_returns_test - current_risk_free_rate;
    
    % Save portfolio return for current week
    portfolio_return_all(week-num_weeks/2) = portfolio_return_test;
    
    % Calculate maximum drawdown for test data using CDD_Var_p function
    max_drawdown_test = CDD_Var_p_cal(week-num_weeks/2, asset_returns, portfolio_weights_all, theta);
    
    % Calculate CDD ratio for test data using portfolio return and maximum drawdown
    CDD_ratio_test = sum(portfolio_return_all) / (week-num_weeks/2) / max_drawdown_test;
    
    % Save CDD ratio for current week
    CDD_ratio_p_all(week-num_weeks/2) = CDD_ratio_test;  
end

% Calculate CDD ratio for test data using different factors depending on the dataset
CDD_ratio_update = CDD_ratio_p_all* (num_weeks/18)^0.5;