% Clear workspace and command window
clc,clear

% Set display format to long
format long;

% Load data from file
direc = './UPDATE_HSCI_Calmar/';
flist = dir([direc, '*.txt']);
asset_prices = load('wk_price_HSCI', ',')';
index_prices = load('b1_HSCI.txt', ','); %index

risk_free_rate=load('rf_00_17.txt',','); %risk-free

% Convert prices to returns
[asset_returns, ~] = price2ret(asset_prices', [], 'Periodic'); %simple return
asset_returns = asset_returns';
[index_returns, ~] = price2ret(index_prices, [], 'Periodic'); %simple return

% Get number of assets and weeks
[num_assets, num_weeks] = size(asset_returns);
risk_free_test = risk_free_rate(num_weeks/2+1:num_weeks);

% Variable initial
Curr_week_temp = 0; 
Actual_week_rt_temp = 0;
Equal_week_rt_temp = 0;
Equal_actual_week_rt_temp = 0;
portfolio_weights_all = zeros(num_assets, num_weeks/2);
equal_weights = 1 / num_assets * ones(num_assets, 1);
equal_weights_all = 1 / num_assets * ones(num_assets, num_weeks/2);
ratio_ew = zeros(1, num_weeks/2);
Curr_week_rt = zeros(1, num_weeks/2);
Actual_week_rt = zeros(1, num_weeks/2);
Equal_week_rt = zeros(1, num_weeks/2);
Equal_week_wk_rt = zeros(1, num_weeks/2);
ratio_p_all = zeros(1, num_weeks/2);
ratio_a_all = zeros(1, num_weeks/2);

for week = (num_weeks/2+1):num_weeks 
    disp(week)
    wk_return_d1_train = asset_returns(:, 1:week-1); %training
    wk_return_d1_test = asset_returns(:, week); 
    rf = risk_free_rate(week); 
    filename = flist(week-num_weeks/2);
    portfolio_weights = load(['./UPDATE_HSCI_Calmar/',filename.name]);
    portfolio_weights_all(:,week-num_weeks/2) = portfolio_weights;
    
    %calmar_ratio_p
    Curr_week_temp = portfolio_weights' * wk_return_d1_test - rf;
    Curr_week_rt(week-num_weeks/2) = Curr_week_temp;
    max_drawdown0 = Calmar_Var_p(week-num_weeks/2, asset_returns, portfolio_weights_all);
    ratio_my = sum(Curr_week_rt) / (week-num_weeks/2) / max_drawdown0;
    ratio_p_all(week-num_weeks/2) = ratio_my;

    %cumulative return
    Actual_week_rt_temp = Actual_week_rt_temp + portfolio_weights' * wk_return_d1_test;
    Actual_week_rt(week-num_weeks/2) = Actual_week_rt_temp;

    %calmar_ratio_a
    My_wk_rt_a_temp = portfolio_weights' * mean(wk_return_d1_train, 2) - rf;
    ratio_a_my = My_wk_rt_a_temp / Calmar_Var_a(week, asset_returns, portfolio_weights);
    ratio_a_all(week-num_weeks/2) = ratio_a_my;

    %equally weighted
    Equal_week_rt_temp = equal_weights' * wk_return_d1_test - rf;
    Equal_week_rt(week-num_weeks/2) = Equal_week_rt_temp;
    ratio_ew = sum(Equal_week_rt) / (week-num_weeks/2) / Calmar_Var_p(week-num_weeks/2, asset_returns, equal_weights_all);
    ratio_ew(week-num_weeks/2) = ratio_ew;

    Equal_actual_week_rt_temp = Equal_actual_week_rt_temp + equal_weights' * wk_return_d1_test;
    Equal_week_wk_rt(week-num_weeks/2) = Equal_actual_week_rt_temp;
    
end

%index
ratio_index = zeros(1, num_weeks/2);
index_wk_rt = 0;
for week = num_weeks/2+1:num_weeks
    index_wk_rt = index_wk_rt + index_returns(:, week) - risk_free_rate(:, week);
    evar_index = Calmar_Var_p_index(week-num_weeks/2, index_returns);
    ratio_index(:, week-num_weeks/2) = index_wk_rt / (week-num_weeks/2) / evar_index;
end

ratio_p_all_yearly = ratio_p_all * (num_weeks/18)^0.5;
ratio_a_all_yearly = ratio_a_all * (num_weeks/18)^0.5;
actual_wk_rt_yearly = Actual_week_rt * (num_weeks/18)^0.5;
ratio_ew_yearly = ratio_ew * (num_weeks/18)^0.5;
ratio_index_yearly = ratio_index * (num_weeks/18)^0.5;


%Display output
ratio_p_all_yearly(end)
ratio_a_all_yearly(end)
actual_wk_rt_yearly(end)
ratio_ew_yearly(end)
ratio_index_yearly(end)


%plotting
time = 1:N/2;
color=[1,0.84314,0;
       0.80392,0.36078,0.36078;
       0,0.749021,1;
       0.82353,0.41176,0.11765;
       0.6902,0.87843,0.90196;
       1,0,0;
       0,1,0;
       0,0,1;];
figure(3);
set(gcf, 'unit', 'centimeters', 'position', [10, 10, 16, 10])
marksize = 5;
linewidth_gca = 1.6;
fontsize_gca = 24;
fontsize_label = 28;
fontsize_legend = 20;
plot(time, ratio_a_all(1:N/2), '-', 'linewidth', 1.6, 'markersize', marksize, 'color', [0 0.80784 0.81961]);
hold on
plot(time, ratio_p_all(1:N/2), '--', 'linewidth', 1.6, 'markersize', marksize, 'color', '#B8860B');
hold on
plot(time, esr_ew(1:N/2), '-.', 'linewidth', 1.6, 'markersize', marksize, 'color', 'r');
hold on
plot(time, esr_index(1:N/2), ':', 'linewidth', 1.6, 'markersize', marksize, 'color', 'm');


xlim([1 N/2])
ylim([0 0.03])
xticks([52 157 261 365 469])
xticklabels({'2010', '2012', '2014', '2016' '2018'})
h = legend('NO$_{\rm \bf a}$', 'NO$_{\rm \bf p}$', 'EW', 'MI');
set(h,'fontsize',fontsize_legend, 'interpreter', 'latex','Location','NorthWest');
set(gca, 'linewidth', linewidth_gca, 'fontsize', fontsize_gca);
ylabel('CR$_{\rm \bf a}$vs CR$_{\rm \bf p}$', 'interpreter', 'latex', 'fontsize', fontsize_label)   ;
title('HSCI'); 
hold on

hfig3 = figure(3);
figWidth = 16;  
figHeight = 10;  
set(hfig3,'PaperUnits','centimeters');
set(hfig3,'PaperPosition',[0 0 figWidth figHeight]);

filefigures = './figures_new/';
saveas(gcf, [filefigures, 'UPDATE_HSCI_CR'], 'eps');
saveas(gcf, [filefigures, 'UPDATE_HSCI_CR'], 'png');