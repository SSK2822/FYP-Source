clc,clear

filedir = './NEW_FTSE100_ESR/';
filelist = dir([filedir, '*.txt']);
rhodir = './NEW_FTSE100_ESR_rho/';
rholist = dir([rhodir, '*.txt']);

load('FTSE100.mat'); %stock price matrix(939 weeks)
wk_return_d1 = Assets_Returns(1:716, :)';
wk_return_b1 = Index_Returns(1:716, :)';
[M, N] = size(wk_return_d1);
rf_02_16=load('rf_02_16.txt',','); %risk free
rt_test = rf_02_16(N/2+1:N);

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

%solutions
disp("Solution started")
%CZeSD
x_czesd = load('./FTSE100/OptPortfolios_CZeSD_FTSE100.txt');
xt_czesd = zeros(M, N);
for ii = 1:N
    if ii>52
        xt_czesd(:, ii) = x_czesd(:, ceil((ii-52)/12));
    end
end

%KP_SSD
x_kpssd = load('./FTSE100/OptPortfolios_KP_SSD_FTSE100.txt');
xt_kpssd = zeros(M, N);
for ii = 1:N
    if ii>52
        xt_kpssd(:, ii) = x_kpssd(:, ceil((ii-52)/12));
    end
end

%L_SSD
x_lssd = load('./FTSE100/OptPortfolios_L_SSD_FTSE100.txt');
xt_lssd = zeros(M, N);
for ii = 1:N
    if ii>52
        xt_lssd(:, ii) = x_lssd(:, ceil((ii-52)/12));
    end
end

%LR_ASSD
x_lrassd = load('./FTSE100/OptPortfolios_LR_ASSD_FTSE100.txt');
xt_lrassd = zeros(M, N);
for ii = 1:N
    if ii>52
        xt_lrassd(:, ii) = x_lrassd(:, ceil((ii-52)/12));
    end
end

%Mean-Variance
x_mv = load('./FTSE100/OptPortfolios_MeanVar_FTSE100.txt');
xt_mv = zeros(M, N);
for ii = 1:N
    if ii>52
        xt_mv(:, ii) = x_mv(:, ceil((ii-52)/12));
    end
end

%RMZ_SSD
x_rmzssd = load('./FTSE100/OptPortfolios_RMZ_SSD_FTSE100.txt');
xt_rmzssd = zeros(M, N);
for ii = 1:N
    if ii>52
        xt_rmzssd(:, ii) = x_rmzssd(:, ceil((ii-52)/12));
    end
end
disp("Solution ended")
disp("NO Started")
for i = (N/2+1):N 
    wk_return_d1_train = wk_return_d1(:, 1:i-1); %training
    wk_return_d1_test = wk_return_d1(:, i); 
    rf = rf_02_16(i); 
    miu = mean(wk_return_d1_train, 2);
    filename = filelist(i-N/2);
    xt = load([filedir,filename.name]);
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
        rhot = load([rhodir,filename_rho.name]);
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
disp("NO Ended")
disp("Index Started")
%index

esr_index = zeros(1, N/2);
index_wk_rt = 0;
for i = N/2+1:N
    index_wk_rt = index_wk_rt + wk_return_b1(:, i) - rf_02_16(:, i);
    [~, evar_index, a] = fmincon(@(rho) EVaR_p_index(i, wk_return_b1, rho, theta), ...
    0.01, [], [], [], [], 0, [], []);
    if (a ~= 1)
        disp(i)
    end
    esr_index(:, i-N/2) = index_wk_rt / (i-N/2) / evar_index;
end
disp("Index Ended")
disp("MV Started")
xt_mv_all = zeros(M, N/2);
esr_mv = zeros(1, N/2);
for ii = (N/2+1):N 
    wk_return_d1_train = wk_return_d1(:, 1:ii-1); %training
    wk_return_d1_test = wk_return_d1(:, ii); 
    rf = rf_02_16(ii); 
    miu = mean(wk_return_d1_train, 2);
    xt_mv_all(:,ii-N/2) = xt_mv(:, ii);
    
    %ESR_p
    My_wk_rt_temp = xt_mv(:, ii)' * wk_return_d1_test - rf;
    My_wk_rt(ii-N/2) = My_wk_rt_temp;
    evar0 = EVaR_p_cal(ii-N/2, wk_return_d1, xt_mv_all, theta);
    ratio_my = sum(My_wk_rt) / (ii-N/2) / evar0;
    esr_mv(ii-N/2) = ratio_my;
end
disp("MV Ended")
disp("LSSD Started")
xt_lssd_all = zeros(M, N/2);
esr_lssd = zeros(1, N/2);
for ii = (N/2+1):N 
    wk_return_d1_train = wk_return_d1(:, 1:ii-1); %training
    wk_return_d1_test = wk_return_d1(:, ii); 
    rf = rf_02_16(ii); 
    miu = mean(wk_return_d1_train, 2);
    xt_lssd_all(:,ii-N/2) = xt_lssd(:, ii);
    
    %ESR_p
    My_wk_rt_temp = xt_lssd(:, ii)' * wk_return_d1_test - rf;
    My_wk_rt(ii-N/2) = My_wk_rt_temp;
    evar0 = EVaR_p_cal(ii-N/2, wk_return_d1, xt_lssd_all, theta);
    ratio_my = sum(My_wk_rt) / (ii-N/2) / evar0;
    esr_lssd(ii-N/2) = ratio_my;
end
disp("LSSD Ended")
disp("LRASSD Started")
xt_lrassd_all = zeros(M, N/2);
esr_lrassd = zeros(1, N/2);
for ii = (N/2+1):N 
    wk_return_d1_train = wk_return_d1(:, 1:ii-1); %training
    wk_return_d1_test = wk_return_d1(:, ii); 
    rf = rf_02_16(ii); 
    miu = mean(wk_return_d1_train, 2);
    xt_lrassd_all(:,ii-N/2) = xt_lrassd(:, ii);
    
    %ESR_p
    My_wk_rt_temp = xt_lrassd(:, ii)' * wk_return_d1_test - rf;
    My_wk_rt(ii-N/2) = My_wk_rt_temp;
    evar0 = EVaR_p_cal(ii-N/2, wk_return_d1, xt_lrassd_all, theta);
    ratio_my = sum(My_wk_rt) / (ii-N/2) / evar0;
    esr_lrassd(ii-N/2) = ratio_my;
end
disp("LRASSD Ended")
disp("RMZSSD Started")
xt_rmzssd_all = zeros(M, N/2);
esr_rmzssd = zeros(1, N/2);
for ii = (N/2+1):N 
    wk_return_d1_train = wk_return_d1(:, 1:ii-1); %training
    wk_return_d1_test = wk_return_d1(:, ii); 
    rf = rf_02_16(ii); 
    miu = mean(wk_return_d1_train, 2);
    xt_rmzssd_all(:,ii-N/2) = xt_rmzssd(:, ii);
    
    %ESR_p
    My_wk_rt_temp = xt_rmzssd(:, ii)' * wk_return_d1_test - rf;
    My_wk_rt(ii-N/2) = My_wk_rt_temp;
    evar0 = EVaR_p_cal(ii-N/2, wk_return_d1, xt_rmzssd_all, theta);
    ratio_my = sum(My_wk_rt) / (ii-N/2) / evar0;
    esr_rmzssd(ii-N/2) = ratio_my;
end
disp("RMZSSD Ended")
disp("KPSSD Started")
xt_kpssd_all = zeros(M, N/2);
esr_kpssd = zeros(1, N/2);
for ii = (N/2+1):N 
    wk_return_d1_train = wk_return_d1(:, 1:ii-1); %training
    wk_return_d1_test = wk_return_d1(:, ii); 
    rf = rf_02_16(ii); 
    miu = mean(wk_return_d1_train, 2);
    xt_kpssd_all(:,ii-N/2) = xt_kpssd(:, ii);
    
    %ESR_p
    My_wk_rt_temp = xt_kpssd(:, ii)' * wk_return_d1_test - rf;
    My_wk_rt(ii-N/2) = My_wk_rt_temp;
    evar0 = EVaR_p_cal(ii-N/2, wk_return_d1, xt_kpssd_all, theta);
    ratio_my = sum(My_wk_rt) / (ii-N/2) / evar0;
    esr_kpssd(ii-N/2) = ratio_my;
end
disp("KPSSD Ended")
disp("CZESD Started")
xt_czesd_all = zeros(M, N/2);
esr_czesd = zeros(1, N/2);
for ii = (N/2+1):N 
    wk_return_d1_train = wk_return_d1(:, 1:ii-1); %training
    wk_return_d1_test = wk_return_d1(:, ii); 
    rf = rf_02_16(ii); 
    miu = mean(wk_return_d1_train, 2);
    xt_czesd_all(:,ii-N/2) = xt_czesd(:, ii);
    
    %ESR_p
    My_wk_rt_temp = xt_czesd(:, ii)' * wk_return_d1_test - rf;
    My_wk_rt(ii-N/2) = My_wk_rt_temp;
    evar0 = EVaR_p_cal(ii-N/2, wk_return_d1, xt_czesd_all, theta);
    ratio_my = sum(My_wk_rt) / (ii-N/2) / evar0;
    esr_czesd(ii-N/2) = ratio_my;
end
disp("CZESD Ended")
esr_yearly = esr * (N/13.83)^0.5;
esr_ew_yearly = esr_ew * (N/13.83)^0.5;

esr_index_yearly = esr_index * (N/13.83)^0.5;
esr_mv_yearly = esr_mv * (N/13.83)^0.5;
esr_lssd_yearly = esr_lssd * (N/13.83)^0.5;
esr_lrassd_yearly = esr_lrassd * (N/13.83)^0.5;
esr_rmzssd_yearly = esr_rmzssd * (N/13.83)^0.5;
esr_kpssd_yearly = esr_kpssd * (N/13.83)^0.5;
esr_czesd_yearly = esr_czesd * (N/13.83)^0.5;

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
ylim([0 0.15])
xticks([28 132 237 341])
xticklabels({'2010', '2012', '2014', '2016'})
h = legend('NO$_{\rm \bf a}$', 'NO$_{\rm \bf p}$', 'EW', 'MI');
set(h,'fontsize',fontsize_legend, 'interpreter', 'latex');
set(gca, 'linewidth', linewidth_gca, 'fontsize', fontsize_gca);
%set(gca, 'GridLineStyle', '--');
ylabel('ESR$_{\rm \bf p}$vs ESR$_{\rm \bf a}$', 'interpreter', 'latex', 'fontsize', fontsize_label)   ;
title('FTSE 100'); % 标题
hold on

hfig3 = figure(3);
figWidth = 16;  % 设置图片宽度
figHeight = 10;  % 设置图片高度
set(hfig3,'PaperUnits','centimeters'); % 图片尺寸所用单位
set(hfig3,'PaperPosition',[0 0 figWidth figHeight]);

%各个esr的最终值比较
disp('solution EW IDX')
esr_yearly(end)
esr_ew_yearly(end)

esr_index_yearly(end)

disp('solution MV LSSD LRASSD RMZSSD KPSSD CZESD')
esr_mv_yearly(end)
esr_lssd_yearly(end)
esr_lrassd_yearly(end)
esr_rmzssd_yearly(end)
esr_kpssd_yearly(end)
esr_czesd_yearly(end)
