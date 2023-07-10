clc,clear

filedir = './NEW_NASDAQ100_Calmar/';
filelist = dir([filedir, '*.txt']);

load('NASDAQ100.mat'); %stock price matrix(939 weeks)
wk_return_d1 = Assets_Returns(1:596, :)';
wk_return_b1 = Index_Returns(1:596, :)';
[M, N] = size(wk_return_d1);
rf_04_16=load('rf_04_16a.txt',','); %risk free
rt_test = rf_04_16(N/2+1:N);

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
ratio_p_all = zeros(1, N/2);
ratio_a_all = zeros(1, N/2);

%solutions
disp("Solution started")
%CZeSD
x_czesd = load('../data/BCST_Solutions/NASDAQ100/OptPortfolios_CZeSD_NASDAQ100.txt');
xt_czesd = zeros(M, N);
for ii = 1:N
    if ii>52
        xt_czesd(:, ii) = x_czesd(:, ceil((ii-52)/12));
    end
end

%KP_SSD
x_kpssd = load('../data/BCST_Solutions/NASDAQ100/OptPortfolios_KP_SSD_NASDAQ100.txt');
xt_kpssd = zeros(M, N);
for ii = 1:N
    if ii>52
        xt_kpssd(:, ii) = x_kpssd(:, ceil((ii-52)/12));
    end
end

%L_SSD
x_lssd = load('../data/BCST_Solutions/NASDAQ100/OptPortfolios_L_SSD_NASDAQ100.txt');
xt_lssd = zeros(M, N);
for ii = 1:N
    if ii>52
        xt_lssd(:, ii) = x_lssd(:, ceil((ii-52)/12));
    end
end

%LR_ASSD
x_lrassd = load('../data/BCST_Solutions/NASDAQ100/OptPortfolios_LR_ASSD_NASDAQ100.txt');
xt_lrassd = zeros(M, N);
for ii = 1:N
    if ii>52
        xt_lrassd(:, ii) = x_lrassd(:, ceil((ii-52)/12));
    end
end

%Mean-Variance
disp("MV Reached")
x_mv = load('../data/BCST_Solutions/NASDAQ100/OptPortfolios_MeanVar_NASDAQ100.txt');
xt_mv = zeros(M, N);
for ii = 1:N
    if ii>52
        xt_mv(:, ii) = x_mv(:, ceil((ii-52)/12));
    end
end

%RMZ_SSD
x_rmzssd = load('../data/BCST_Solutions/NASDAQ100/OptPortfolios_RMZ_SSD_NASDAQ100.txt');
xt_rmzssd = zeros(M, N);
for ii = 1:N
    if ii>52
        xt_rmzssd(:, ii) = x_rmzssd(:, ceil((ii-52)/12));
    end
end
disp("Solution ended")
for i = (N/2+1):N 
    i
    wk_return_d1_train = wk_return_d1(:, 1:i-1); %training
    wk_return_d1_test = wk_return_d1(:, i); 
    rf = rf_04_16(i); 
    miu = mean(wk_return_d1_train, 2);
    filename = filelist(i-N/2);
    xt = load([filedir,filename.name]);
    xt_all(:,i-N/2) = xt;
    
    %ESR_p
    My_wk_rt_temp = xt' * wk_return_d1_test - rf;
    My_wk_rt(i-N/2) = My_wk_rt_temp;
    evar0 = Calmar_Var_p(i-N/2, wk_return_d1, xt_all);
    ratio_my = sum(My_wk_rt) / (i-N/2) / evar0;
    esr(i-N/2) = ratio_my;

    %cumulative return
    ac_wk_rt_temp = ac_wk_rt_temp + xt' * wk_return_d1_test;
    ac_wk_rt(i-N/2) = ac_wk_rt_temp;

    %esr_a
    
    My_wk_rt_a_temp = xt' * mean(wk_return_d1_train, 2) - rf;
    ratio_a_my = My_wk_rt_a_temp / Calmar_Var_a(i, wk_return_d1, xt);
    esr_a(i-N/2) = ratio_a_my;
  

    %equally weighted
    ew_wk_rt_temp = x_ew' * wk_return_d1_test - rf;
    ew_wk_rt(i-N/2) = ew_wk_rt_temp;
    ratio_ew = sum(ew_wk_rt) / (i-N/2) / Calmar_Var_p(i-N/2, wk_return_d1, x_ew_all);
    esr_ew(i-N/2) = ratio_ew;
    
    ew_ac_wk_rt_temp = ew_ac_wk_rt_temp + x_ew' * wk_return_d1_test;
    ew_ac_wk_rt(i-N/2) = ew_ac_wk_rt_temp;
    
end

%index

esr_index = zeros(1, N/2);%初始化
index_wk_rt = 0;
for i = N/2+1:N
    index_wk_rt = index_wk_rt + wk_return_b1(:, i) - rf_04_16(:, i);
    %[~, evar_index, a] = fmincon(@(rho) EVaR_p_index(i, wk_return_b1, rho, theta), ...
    %0.01, [], [], [], [], 0, [], []);
    %if (a ~= 1)
     %   disp(i)
    %end
    evar_index = Calmar_Var_p_index(i-N/2, wk_return_b1);
    esr_index(:, i-N/2) = index_wk_rt / (i-N/2) / evar_index;
end



xt_mv_all = zeros(M, N/2);
esr_mv = zeros(1, N/2);
for ii = (N/2+1):N 
    wk_return_d1_train = wk_return_d1(:, 1:ii-1); %training
    wk_return_d1_test = wk_return_d1(:, ii); 
    rf = rf_04_16(ii); 
    miu = mean(wk_return_d1_train, 2);
    xt_mv_all(:,ii-N/2) = xt_mv(:, ii);
    
    %ESR_p
    My_wk_rt_temp = xt_mv(:, ii)' * wk_return_d1_test - rf;
    My_wk_rt(ii-N/2) = My_wk_rt_temp;
    evar0 = Calmar_Var_p(ii-N/2, wk_return_d1, xt_mv_all);
    ratio_my = sum(My_wk_rt) / (ii-N/2) / evar0;
    esr_mv(ii-N/2) = ratio_my;
end

xt_lssd_all = zeros(M, N/2);
esr_lssd = zeros(1, N/2);
for ii = (N/2+1):N 
    wk_return_d1_train = wk_return_d1(:, 1:ii-1); %training
    wk_return_d1_test = wk_return_d1(:, ii); 
    rf = rf_04_16(ii); 
    miu = mean(wk_return_d1_train, 2);
    xt_lssd_all(:,ii-N/2) = xt_lssd(:, ii);
    
    %ESR_p
    My_wk_rt_temp = xt_lssd(:, ii)' * wk_return_d1_test - rf;
    My_wk_rt(ii-N/2) = My_wk_rt_temp;
    evar0 = Calmar_Var_p(ii-N/2, wk_return_d1, xt_lssd_all);
    ratio_my = sum(My_wk_rt) / (ii-N/2) / evar0;
    esr_lssd(ii-N/2) = ratio_my;
end

xt_lrassd_all = zeros(M, N/2);
esr_lrassd = zeros(1, N/2);
for ii = (N/2+1):N 
    wk_return_d1_train = wk_return_d1(:, 1:ii-1); %training
    wk_return_d1_test = wk_return_d1(:, ii); 
    rf = rf_04_16(ii); 
    miu = mean(wk_return_d1_train, 2);
    xt_lrassd_all(:,ii-N/2) = xt_lrassd(:, ii);
    
    %ESR_p
    My_wk_rt_temp = xt_lrassd(:, ii)' * wk_return_d1_test - rf;
    My_wk_rt(ii-N/2) = My_wk_rt_temp;
    evar0 = Calmar_Var_p(ii-N/2, wk_return_d1, xt_lrassd_all);
    ratio_my = sum(My_wk_rt) / (ii-N/2) / evar0;
    esr_lrassd(ii-N/2) = ratio_my;
end

xt_rmzssd_all = zeros(M, N/2);
esr_rmzssd = zeros(1, N/2);
for ii = (N/2+1):N 
    wk_return_d1_train = wk_return_d1(:, 1:ii-1); %training
    wk_return_d1_test = wk_return_d1(:, ii); 
    rf = rf_04_16(ii); 
    miu = mean(wk_return_d1_train, 2);
    xt_rmzssd_all(:,ii-N/2) = xt_rmzssd(:, ii);
    
    %ESR_p
    My_wk_rt_temp = xt_rmzssd(:, ii)' * wk_return_d1_test - rf;
    My_wk_rt(ii-N/2) = My_wk_rt_temp;
    evar0 = Calmar_Var_p(ii-N/2, wk_return_d1, xt_rmzssd_all);
    ratio_my = sum(My_wk_rt) / (ii-N/2) / evar0;
    esr_rmzssd(ii-N/2) = ratio_my;
end

xt_kpssd_all = zeros(M, N/2);
esr_kpssd = zeros(1, N/2);
for ii = (N/2+1):N 
    wk_return_d1_train = wk_return_d1(:, 1:ii-1); %training
    wk_return_d1_test = wk_return_d1(:, ii); 
    rf = rf_04_16(ii); 
    miu = mean(wk_return_d1_train, 2);
    xt_kpssd_all(:,ii-N/2) = xt_kpssd(:, ii);
    
    %ESR_p
    My_wk_rt_temp = xt_kpssd(:, ii)' * wk_return_d1_test - rf;
    My_wk_rt(ii-N/2) = My_wk_rt_temp;
    evar0 = Calmar_Var_p(ii-N/2, wk_return_d1, xt_kpssd_all);
    ratio_my = sum(My_wk_rt) / (ii-N/2) / evar0;
    esr_kpssd(ii-N/2) = ratio_my;
end

xt_czesd_all = zeros(M, N/2);
esr_czesd = zeros(1, N/2);
for ii = (N/2+1):N 
    wk_return_d1_train = wk_return_d1(:, 1:ii-1); %training
    wk_return_d1_test = wk_return_d1(:, ii); 
    rf = rf_04_16(ii); 
    miu = mean(wk_return_d1_train, 2);
    xt_czesd_all(:,ii-N/2) = xt_czesd(:, ii);
    
    %ESR_p
    My_wk_rt_temp = xt_czesd(:, ii)' * wk_return_d1_test - rf;
    My_wk_rt(ii-N/2) = My_wk_rt_temp;
    evar0 = Calmar_Var_p(ii-N/2, wk_return_d1, xt_czesd_all);
    ratio_my = sum(My_wk_rt) / (ii-N/2) / evar0;
    esr_czesd(ii-N/2) = ratio_my;
end

esr_yearly = esr * (N/11.5)^0.5;
esr_ew_yearly = esr_ew * (N/11.5)^0.5;

esr_index_yearly = esr_index * (N/11.5)^0.5;
esr_mv_yearly = esr_mv * (N/11.5)^0.5;
esr_lssd_yearly = esr_lssd * (N/11.5)^0.5;
esr_lrassd_yearly = esr_lrassd * (N/11.5)^0.5;
esr_rmzssd_yearly = esr_rmzssd * (N/11.5)^0.5;
esr_kpssd_yearly = esr_kpssd * (N/11.5)^0.5;
esr_czesd_yearly = esr_czesd * (N/11.5)^0.5;

%开始画图

time = [1:N/2];
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
plot(time, esr_a(1:N/2), '-', 'linewidth', 1.6, 'markersize', marksize, 'color', [0 0.80784 0.81961]);
hold on
plot(time, esr(1:N/2), '--', 'linewidth', 1.6, 'markersize', marksize, 'color', '#B8860B');
hold on
plot(time, esr_ew(1:N/2), '-.', 'linewidth', 1.6, 'markersize', marksize, 'color', 'r');
hold on
plot(time, esr_index(1:N/2), ':', 'linewidth', 1.6, 'markersize', marksize, 'color', 'm');
xlim([1 N/2])
ylim([0 0.14])
xticks([20 72 124 177 229 281])
xticklabels({'2011', '2012', '2013', '2014', '2015', '2016'})
h = legend('NO$_{\rm \bf a}$', 'NO$_{\rm \bf p}$', 'EW', 'MI');
set(h,'fontsize',fontsize_legend, 'interpreter', 'latex');
set(gca, 'linewidth', linewidth_gca, 'fontsize', fontsize_gca);
%set(gca, 'GridLineStyle', '--');
ylabel('ESR$_{\rm \bf p}$vs ESR$_{\rm \bf a}$', 'interpreter', 'latex', 'fontsize', fontsize_label)   ;
title('NASDAQ 100'); % 标题
hold on

hfig3 = figure(3);
figWidth = 16;  % 设置图片宽度
figHeight = 10;  % 设置图片高度
set(hfig3,'PaperUnits','centimeters'); % 图片尺寸所用单位
set(hfig3,'PaperPosition',[0 0 figWidth figHeight]);

%各个esr的最终值比较
esr_yearly(end)
esr_ew_yearly(end)

esr_index_yearly(end)

disp('solution')
esr_mv_yearly(end)
esr_lssd_yearly(end)
esr_lrassd_yearly(end)
esr_rmzssd_yearly(end)
esr_kpssd_yearly(end)
esr_czesd_yearly(end)
