clc,clear

filedir = './UPDATE_SP500_ESR/';
filelist = dir([filedir, '*.txt']);
rhodir = './UPDATE_SP500_ESR_rho/';
rholist = dir([rhodir, '*.txt']);
d1 = load('wk_price_SP500', ',')';
b1 = load('SP500_historical_index', ',')'; %index
rf_00_17=load('rf_00_17.txt',','); %risk-free


[wk_return_d1, ~] = price2ret(d1', [], 'Periodic'); %simple return
wk_return_d1 = wk_return_d1';
[wk_return_b1, ~] = price2ret(b1, [], 'Periodic'); %simple return

[M, N] = size(wk_return_d1);
rf_test = rf_00_17(N/2+1:N);

My_wk_rt_temp = 0; 
ac_wk_rt_temp = 0;
ew_wk_rt_temp = 0;
ew_ac_wk_rt_temp = 0;
xt_all = zeros(M, N/2);
x_ew = 1 / M * ones(M, 1);
x_ew_all = 1 / M * ones(M, N/2);
esr = zeros(1, N/2);
esr_ew = zeros(1, N/2);
My_wk_rt = zeros(1, N/2);
ac_wk_rt = zeros(1, N/2);
ew_wk_rt = zeros(1, N/2);
ew_ac_wk_rt = zeros(1, N/2);
esr_a = zeros(1, N/2);
theta = 0.95;

for i = (N/2+1):N 
    wk_return_d1_train = wk_return_d1(:, 1:i-1); %training
    wk_return_d1_test = wk_return_d1(:, i); 
    rf = rf_00_17(i); 
    miu = mean(wk_return_d1_train, 2);
    filename = filelist(i-N/2);
    xt = load(['./UPDATE_SP500_ESR/',filename.name]);
    xt_all(:,i-N/2) = xt;
    
    %ESR_p
    My_wk_rt_temp = xt' * wk_return_d1_test - rf;
    My_wk_rt(i-N/2) = My_wk_rt_temp;
    evar0 = EVaR_p_cal(i-N/2, wk_return_d1, xt_all, theta);
    ratio_my = sum(My_wk_rt) / (i-N/2) / evar0;
    esr(i-N/2) = ratio_my;

    %cumulative return
    ac_wk_rt_temp = ac_wk_rt_temp + xt' * wk_return_d1_test;
    ac_wk_rt(i-N/2) = ac_wk_rt_temp;

    %esr_a
    if (i<N)
        filename_rho = rholist(i-N/2+1);
        rhot = load(['./UPDATE_SP500_ESR_rho/',filename_rho.name]);
        My_wk_rt_a_temp = xt' * mean(wk_return_d1_train, 2) - rf;
        ratio_a_my = My_wk_rt_a_temp / EVaR_a(i, wk_return_d1, xt, rhot, theta);
        esr_a(i-N/2) = ratio_a_my;
    end

    if (i == N)
        My_wk_rt_a_temp = xt' * mean(wk_return_d1_train, 2) - rf;
        ratio_a_my = My_wk_rt_a_temp / EVaR_a_cal(i, wk_return_d1, xt, theta);
        esr_a(i-N/2) = ratio_a_my;
    end

    %equally weighted
    ew_wk_rt_temp = x_ew' * wk_return_d1_test - rf;
    ew_wk_rt(i-N/2) = ew_wk_rt_temp;
    ratio_ew = sum(ew_wk_rt) / (i-N/2) / EVaR_p_cal(i-N/2, wk_return_d1, x_ew_all, theta);
    esr_ew(i-N/2) = ratio_ew;
    
    ew_ac_wk_rt_temp = ew_ac_wk_rt_temp + x_ew' * wk_return_d1_test;
    ew_ac_wk_rt(i-N/2) = ew_ac_wk_rt_temp;
    
end

%index
esr_index = zeros(1, N/2);%初始化
index_wk_rt = 0;
for i = N/2+1:N
    index_wk_rt = index_wk_rt + wk_return_b1(:, i) - rf_00_17(:, i);
    [~, evar_index, a] = fmincon(@(rho) EVaR_p_index(i, wk_return_b1, rho, theta), ...
    10, [], [], [], [], 0, [], []);
    disp(['index_p_cal:' num2str(a)]);
    esr_index(:, i-N/2) = index_wk_rt / (i-N/2) / evar_index;
end

esr_yearly = esr * (938/18)^0.5;
esr_ew_yearly = esr_ew * (938/18)^0.5;
esr_index_yearly = esr_index * (938/18)^0.5;


%开始画图
%{
time = [1:469];
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
plot(time, esr_a(1:469), '-', 'linewidth', 1.6, 'markersize', marksize, 'color', [0 0.80784 0.81961]);
hold on
plot(time, esr(1:469), '--', 'linewidth', 1.6, 'markersize', marksize, 'color', '#B8860B');
hold on
plot(time, esr_ew(1:469), '-.', 'linewidth', 1.6, 'markersize', marksize, 'color', 'r');
hold on
plot(time, esr_index(1:469), ':', 'linewidth', 1.6, 'markersize', marksize, 'color', 'm');
xlim([1 469])
ylim([0 0.2])
xticks([52 157 261 365 469])
xticklabels({'2010', '2012', '2014', '2016', '2018'})
h = legend('NO$_{\rm \bf a}$', 'NO$_{\rm \bf p}$', 'EW', 'MI');
set(h,'fontsize',fontsize_legend, 'interpreter', 'latex');
set(gca, 'linewidth', linewidth_gca, 'fontsize', fontsize_gca);
%set(gca, 'GridLineStyle', '--');
ylabel('ESR$_{\rm \bf p}$vs ESR$_{\rm \bf a}$', 'interpreter', 'latex', 'fontsize', fontsize_label)   ;
title('SP 500'); % 标题
hold on

hfig3 = figure(3);
figWidth = 16;  % 设置图片宽度
figHeight = 10;  % 设置图片高度
set(hfig3,'PaperUnits','centimeters'); % 图片尺寸所用单位
set(hfig3,'PaperPosition',[0 0 figWidth figHeight]);
%}

%各个esr的最终值比较
esr_yearly(end)
esr_ew_yearly(end)
esr_index_yearly(end)