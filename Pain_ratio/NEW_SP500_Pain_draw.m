% Clear workspace and command window
clc,clear

% Set display format to long
format long;
filedir = './NEW_SP500_Pain/';
filelist = dir([filedir, '*.txt']);

% Load data from file
load('SP500.mat'); %stock price matrix(939 weeks)
asset_returns = Assets_Returns(1:1362, :)';
index_returns = Index_Returns(1:1362, :)';

[num_assets, num_weeks] = size(asset_returns);
risk_free_rate=load('rf_04_16b.txt',','); %risk free
rt_test = risk_free_rate(num_weeks/2+1:num_weeks);

Curr_week_temp = 0; 
Actual_week_rt_temp = 0;
Equal_week_rt_temp = 0;
Equal_actual_week_rt_temp = 0;
portfolio_weights_all = zeros(num_assets, num_weeks/2);
equal_weights = 1 / num_assets * ones(num_assets, 1);
equal_weights_all = 1 / num_assets * ones(num_assets, num_weeks/2);
ratio = zeros(1, num_weeks/2);
ratio_ew = zeros(1, num_weeks/2);
Curr_week_rt = zeros(1, num_weeks/2);
Actual_week_rt = zeros(1, num_weeks/2);
Equal_week_rt = zeros(1, num_weeks/2);
Equal_week_wk_rt = zeros(1, num_weeks/2);
ratio_a = zeros(1, num_weeks/2);
ratio_p_all = zeros(1, num_weeks/2);
ratio_a_all = zeros(1, num_weeks/2);

%CZeSD
x_czesd = load('../data/BCST_Solutions/SP500/OptPortfolios_CZeSD_SP500.txt');
portfolio_czesd = zeros(num_assets, num_weeks);
for weeks = 1:num_weeks
    if weeks>52
        portfolio_czesd(:, weeks) = x_czesd(:, ceil((weeks-52)/12));
    end
end

%KP_SSD
x_kpssd = load('../data/BCST_Solutions/SP500/OptPortfolios_KP_SSD_SP500.txt');
portfolio_kpssd = zeros(num_assets, num_weeks);
for weeks = 1:num_weeks
    if weeks>52
        portfolio_kpssd(:, weeks) = x_kpssd(:, ceil((weeks-52)/12));
    end
end

%L_SSD
x_lssd = load('../data/BCST_Solutions/SP500/OptPortfolios_L_SSD_SP500.txt');
portfolio_lssd = zeros(num_assets, num_weeks);
for weeks = 1:num_weeks
    if weeks>52
        portfolio_lssd(:, weeks) = x_lssd(:, ceil((weeks-52)/12));
    end
end

%LR_ASSD
x_lrassd = load('../data/BCST_Solutions/SP500/OptPortfolios_LR_ASSD_SP500.txt');
portfolio_lrassd = zeros(num_assets, num_weeks);
for weeks = 1:num_weeks
    if weeks>52
        portfolio_lrassd(:, weeks) = x_lrassd(:, ceil((weeks-52)/12));
    end
end

%MV
x_mv = load('../data/BCST_Solutions/SP500/OptPortfolios_MeanVar_SP500.txt');
portfolio_mv = zeros(num_assets, num_weeks);
for weeks = 1:num_weeks
    if weeks>52
        portfolio_mv(:, weeks) = x_mv(:, ceil((weeks-52)/12));
    end
end

%RMZ_SSD
x_rmzssd = load('../data/BCST_Solutions/SP500/OptPortfolios_RMZ_SSD_SP500.txt');
portfolio_rmzssd = zeros(num_assets, num_weeks);
for weeks = 1:num_weeks
    if weeks>52
        portfolio_rmzssd(:, weeks) = x_rmzssd(:, ceil((weeks-52)/12));
    end
end

disp("NO Started")
for week = (num_weeks/2+1):num_weeks 
    disp(week)
    asset_return_d1_train = asset_returns(:, 1:week-1); %training
    asset_return_d1_test = asset_returns(:, week); 
    rf = risk_free_rate(week); 
    miu = mean(asset_return_d1_train, 2);
    filename = filelist(week-num_weeks/2);
    portfolio = load([filedir,filename.name]);
    portfolio_weights_all(:,week-num_weeks/2) = portfolio;
    
    %Pain_ratio_p
    Curr_week_temp = portfolio' * asset_return_d1_test - rf;
    Curr_week_rt(week-num_weeks/2) = Curr_week_temp;
    max_drawdown0 = Pain_Var_p(week-num_weeks/2, asset_returns, portfolio_weights_all);
    ratio_my = sum(Curr_week_rt) / (week-num_weeks/2) / max_drawdown0;
    ratio(week-num_weeks/2) = ratio_my;

    %cumulative return
    Actual_week_rt_temp = Actual_week_rt_temp + portfolio' * asset_return_d1_test;
    Actual_week_rt(week-num_weeks/2) = Actual_week_rt_temp;

    %Pain_ratio_a
    My_wk_rt_a_temp = portfolio' * mean(asset_return_d1_train, 2) - rf;
    ratio_a_my = My_wk_rt_a_temp / Pain_Var_a(week, asset_returns, portfolio);
    ratio_a(week-num_weeks/2) = ratio_a_my;

    %equally weighted
    Equal_week_rt_temp = equal_weights' * asset_return_d1_test - rf;
    Equal_week_rt(week-num_weeks/2) = Equal_week_rt_temp;
    ratio_ew = sum(Equal_week_rt) / (week-num_weeks/2) / Pain_Var_p(week-num_weeks/2, asset_returns, equal_weights_all);
    ratio_ew(week-num_weeks/2) = ratio_ew;
    
    Equal_actual_week_rt_temp = Equal_actual_week_rt_temp + equal_weights' * asset_return_d1_test;
    Equal_week_wk_rt(week-num_weeks/2) = Equal_actual_week_rt_temp;
    
end
disp("NO Ended")
disp("Index Started")
ratio_index = zeros(1, num_weeks/2);
index_wk_rt = 0;
for week = num_weeks/2+1:num_weeks
    index_wk_rt = index_wk_rt + index_returns(:, week) - risk_free_rate(:, week);
    evar_index = Pain_Var_p_index(week-num_weeks/2, index_returns);
    ratio_index(:, week-num_weeks/2) = index_wk_rt / (week-num_weeks/2) / evar_index;
end
disp("Index Ended")
disp("MV Started")
portfolio_mv_all = zeros(num_assets, num_weeks/2);
ratio_mv = zeros(1, num_weeks/2);
for weeks = (num_weeks/2+1):num_weeks 
    disp(weeks)
    asset_return_d1_train = asset_returns(:, 1:weeks-1); %training
    asset_return_d1_test = asset_returns(:, weeks); 
    rf = risk_free_rate(weeks); 
    miu = mean(asset_return_d1_train, 2);
    portfolio_mv_all(:,weeks-num_weeks/2) = portfolio_mv(:, weeks);
    
    Curr_week_temp = portfolio_mv(:, weeks)' * asset_return_d1_test - rf;
    Curr_week_rt(weeks-num_weeks/2) = Curr_week_temp;
    ratio0 = Pain_Var_p(weeks-num_weeks/2, asset_returns, portfolio_mv_all);
    ratio_my = sum(Curr_week_rt) / (weeks-num_weeks/2) / ratio0;
    ratio_mv(weeks-num_weeks/2) = ratio_my;
end
disp("MV Ended")
disp("LSSD Started")
portfolio_lssd_all = zeros(num_assets, num_weeks/2);
ratio_lssd = zeros(1, num_weeks/2);
for weeks = (num_weeks/2+1):num_weeks 
    disp(weeks)
    asset_return_d1_train = asset_returns(:, 1:weeks-1); %training
    asset_return_d1_test = asset_returns(:, weeks); 
    rf = risk_free_rate(weeks); 
    miu = mean(asset_return_d1_train, 2);
    portfolio_lssd_all(:,weeks-num_weeks/2) = portfolio_lssd(:, weeks);

    Curr_week_temp = portfolio_lssd(:, weeks)' * asset_return_d1_test - rf;
    Curr_week_rt(weeks-num_weeks/2) = Curr_week_temp;
    ratio0 = Pain_Var_p(weeks-num_weeks/2, asset_returns, portfolio_lssd_all);
    ratio_my = sum(Curr_week_rt) / (weeks-num_weeks/2) / ratio0;
    ratio_lssd(weeks-num_weeks/2) = ratio_my;
end
disp("LSSD Ended")
disp("LRASSD Started")
portfolio_lrassd_all = zeros(num_assets, num_weeks/2);
ratio_lrassd = zeros(1, num_weeks/2);
for weeks = (num_weeks/2+1):num_weeks 
    disp(weeks)
    asset_return_d1_train = asset_returns(:, 1:weeks-1); %training
    asset_return_d1_test = asset_returns(:, weeks); 
    rf = risk_free_rate(weeks); 
    miu = mean(asset_return_d1_train, 2);
    portfolio_lrassd_all(:,weeks-num_weeks/2) = portfolio_lrassd(:, weeks);
    
    Curr_week_temp = portfolio_lrassd(:, weeks)' * asset_return_d1_test - rf;
    Curr_week_rt(weeks-num_weeks/2) = Curr_week_temp;
    ratio0 = Pain_Var_p(weeks-num_weeks/2, asset_returns, portfolio_lrassd_all);
    ratio_my = sum(Curr_week_rt) / (weeks-num_weeks/2) / ratio0;
    ratio_lrassd(weeks-num_weeks/2) = ratio_my;
end
disp("LRASSD Ended")
disp("RMZSSD Started")
portfolio_rmzssd_all = zeros(num_assets, num_weeks/2);
ratio_rmzssd = zeros(1, num_weeks/2);
for weeks = (num_weeks/2+1):num_weeks 
    disp(weeks)
    asset_return_d1_train = asset_returns(:, 1:weeks-1); %training
    asset_return_d1_test = asset_returns(:, weeks); 
    rf = risk_free_rate(weeks); 
    miu = mean(asset_return_d1_train, 2);
    portfolio_rmzssd_all(:,weeks-num_weeks/2) = portfolio_rmzssd(:, weeks);

    Curr_week_temp = portfolio_rmzssd(:, weeks)' * asset_return_d1_test - rf;
    Curr_week_rt(weeks-num_weeks/2) = Curr_week_temp;
    ratio0 = Pain_Var_p(weeks-num_weeks/2, asset_returns, portfolio_rmzssd_all);
    ratio_my = sum(Curr_week_rt) / (weeks-num_weeks/2) / ratio0;
    ratio_rmzssd(weeks-num_weeks/2) = ratio_my;
end
disp("RMZSSD Ended")
disp("KPSSD Started")
portfolio_kpssd_all = zeros(num_assets, num_weeks/2);
ratio_kpssd = zeros(1, num_weeks/2);
for weeks = (num_weeks/2+1):num_weeks 
    disp(weeks)
    asset_return_d1_train = asset_returns(:, 1:weeks-1); %training
    asset_return_d1_test = asset_returns(:, weeks); 
    rf = risk_free_rate(weeks); 
    miu = mean(asset_return_d1_train, 2);
    portfolio_kpssd_all(:,weeks-num_weeks/2) = portfolio_kpssd(:, weeks);

    Curr_week_temp = portfolio_kpssd(:, weeks)' * asset_return_d1_test - rf;
    Curr_week_rt(weeks-num_weeks/2) = Curr_week_temp;
    ratio0 = Pain_Var_p(weeks-num_weeks/2, asset_returns, portfolio_kpssd_all);
    ratio_my = sum(Curr_week_rt) / (weeks-num_weeks/2) / ratio0;
    ratio_kpssd(weeks-num_weeks/2) = ratio_my;
end
disp("KPSSD Ended")
disp("CZESD Started")
portfolio_czesd_all = zeros(num_assets, num_weeks/2);
ratio_czesd = zeros(1, num_weeks/2);
for weeks = (num_weeks/2+1):num_weeks 
    disp(weeks)
    asset_return_d1_train = asset_returns(:, 1:weeks-1); %training
    asset_return_d1_test = asset_returns(:, weeks); 
    rf = risk_free_rate(weeks); 
    miu = mean(asset_return_d1_train, 2);
    portfolio_czesd_all(:,weeks-num_weeks/2) = portfolio_czesd(:, weeks);

    Curr_week_temp = portfolio_czesd(:, weeks)' * asset_return_d1_test - rf;
    Curr_week_rt(weeks-num_weeks/2) = Curr_week_temp;
    ratio0 = Pain_Var_p(weeks-num_weeks/2, asset_returns, portfolio_czesd_all);
    ratio_my = sum(Curr_week_rt) / (weeks-num_weeks/2) / ratio0;
    ratio_czesd(weeks-num_weeks/2) = ratio_my;
end
disp("CZESD Ended")

ratio_yearly = ratio * (num_weeks/13.83)^0.5;
ratio_ew_yearly = ratio_ew * (num_weeks/13.83)^0.5;
ratio_index_yearly = ratio_index * (num_weeks/13.83)^0.5;
ratio_mv_yearly = ratio_mv * (num_weeks/13.83)^0.5;
ratio_lssd_yearly = ratio_lssd * (num_weeks/13.83)^0.5;
ratio_lrassd_yearly = ratio_lrassd * (num_weeks/13.83)^0.5;
ratio_rmzssd_yearly = ratio_rmzssd * (num_weeks/13.83)^0.5;
ratio_kpssd_yearly = ratio_kpssd * (num_weeks/13.83)^0.5;
ratio_czesd_yearly = ratio_czesd * (num_weeks/13.83)^0.5;

disp('solution EW IDX')
ratio_yearly(end)
ratio_ew_yearly(end)
ratio_index_yearly(end)

disp('solution MV LSSD LRASSD RMZSSD KPSSD CZESD')
ratio_mv_yearly(end)
ratio_lssd_yearly(end)
ratio_lrassd_yearly(end)
ratio_rmzssd_yearly(end)
ratio_kpssd_yearly(end)
ratio_czesd_yearly(end)

%Plotting
time = 1:num_weeks/2;
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
plot(time, ratio_a(1:num_weeks/2), '-', 'linewidth', 1.6, 'markersize', marksize, 'color', [0 0.80784 0.81961]);
hold on
plot(time, ratio(1:num_weeks/2), '--', 'linewidth', 1.6, 'markersize', marksize, 'color', '#B8860B');
hold on
plot(time, ratio_ew(1:num_weeks/2), '-.', 'linewidth', 1.6, 'markersize', marksize, 'color', 'r');
hold on
plot(time, ratio_index(1:num_weeks/2), ':', 'linewidth', 1.6, 'markersize', marksize, 'color', 'm');

xlim([1 N/2])
ylim([0 0.15])
xticks([20 72 124 177 229 281])
xticklabels({'2011', '2012', '2013', '2014', '2015', '2016'})
h = legend('NO$_{\rm \bf a}$', 'NO$_{\rm \bf p}$', 'EW', 'MI');
set(h,'fontsize',fontsize_legend, 'interpreter', 'latex');
set(gca, 'linewidth', linewidth_gca, 'fontsize', fontsize_gca);
ylabel('PR$_{\rm \bf p}$vs PR$_{\rm \bf a}$', 'interpreter', 'latex', 'fontsize', fontsize_label)   ;
title('SP 500');
hold on

hfig3 = figure(3);
figWidth = 16;
figHeight = 10;
set(hfig3,'PaperUnits','centimeters');
set(hfig3,'PaperPosition',[0 0 figWidth figHeight]);