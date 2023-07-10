clc,clear
format long;

filedir = './OR_FTSE100_SOR/';
filelist = dir([filedir, '*.txt']);
d1 = load('ORalldata_FTSE.txt', ','); %price矩阵(939天数据)
b1 = load('ORb1_FTSE.txt', ','); %股指 
rf_92_97 = load('rf_92_97.txt',','); %risk-free

rt_test = rf_92_97(146:290);%后半期risk-free

[wk_return_d1, ~] = price2ret(d1', [], 'Periodic'); %simple return
wk_return_d1 = wk_return_d1';
[wk_return_b1, ~] = price2ret(b1, [], 'Periodic'); %simple return

wk_return_d1 = wk_return_d1(:, 1:290); %weekly return
wk_return_b1 = wk_return_b1(:, 1:290); %指数化为return

[M, N] = size(wk_return_d1);
My_wk_rt_temp = 0; %初始化用于后面ratio分子的计算
ac_wk_rt_temp = 0;
ew_wk_rt_temp = 0;
ew_ac_wk_rt_temp = 0;
xt = zeros(M, 1);%初始化输出列向量xt
xt_all = zeros(M, N/2);
x_ew = 1 / M * ones(M, 1);
sor = zeros(1, 145);%初始化行向量%目标输出测试集的所有sr
sor_ew = zeros(1, 145);
My_wk_rt = zeros(1, 145);
ac_wk_rt = zeros(1, 145);
ew_wk_rt = zeros(1, 145);
ew_ac_wk_rt = zeros(1, 145);
sor_a_od = zeros(1, 145);

for i = (N/2+1):N %后半期的第i周输出：xt，sor_p,绘图：
    wk_return_d1_train = wk_return_d1(:, 1:i-1); %用以训练的return
    wk_return_d1_test = wk_return_d1(:, i); %第i时期的return
    rf = rf_92_97(i); %第i时期的risk free
    miu = mean(wk_return_d1_train, 2); %每个时期的miu %每期更新
    filename = filelist(i-145);
    xt = load(['./OR_FTSE100_SOR/',filename.name]);
    xt_all(:, i-N/2) = xt;
    
    
    %计算sor值  %可以顺便只求return
    My_wk_rt_temp = xt' * wk_return_d1_test - rf;
    My_wk_rt(i-145) = My_wk_rt_temp;
    ac_wk_rt_temp = ac_wk_rt_temp + xt' * wk_return_d1_test;
    ac_wk_rt(i-145) = ac_wk_rt_temp;
    dsr0 = DSR_p(i, wk_return_d1, xt_all, rf_92_97);
    ratio_my = sum(My_wk_rt) / (i-145) / (dsr0 ^ 0.5);
    sor(i-145) = ratio_my;

    %计算先验sor值
    My_wk_rt_a_temp = xt' * mean(wk_return_d1_train, 2) - rf;
    ratio_a_my = My_wk_rt_a_temp / (DSR_a(i, wk_return_d1, xt, rf_92_97))^0.5;
    sor_a_od(i-145) = ratio_a_my;

    %equally weighted sr_p
    ew_wk_rt_temp = x_ew' * wk_return_d1_test - rf;
    ew_wk_rt(i-145) = ew_wk_rt_temp;
    ew_ac_wk_rt_temp = ew_ac_wk_rt_temp + x_ew' * wk_return_d1_test;
    ew_ac_wk_rt(i-145) = ew_ac_wk_rt_temp;
    dsr0 = DSR_p_or(i-145, wk_return_d1, x_ew, rf_92_97);
    ratio_ew = sum(ew_wk_rt) / (i-145) / (dsr0 ^ 0.5);
    sor_ew(i-145) = ratio_ew;
end
%任务1:sor的比较
%index的sor
wk_return_b1_p = wk_return_b1((N/2+1):N);
dsr_index = zeros(1, N/2); %初始化
dsr = zeros(1, N/2);%初始化
dsr_temp = zeros(1, N/2);
for i = 1:(N/2)
    for j = 1:i
        if (wk_return_b1_p(j) < rt_test(j))
            dsr(j) = (wk_return_b1_p(j) - rt_test(j))^2;
            dsr_temp(i) = sum(dsr(1:i));
        end
        dsr_index(i) = dsr_temp(i) / i;
    end
end
sor_index = zeros(1, N/2);
index_ac_wk_rt = zeros(1, N/2);
for i = 1:(N/2)
    wk_return_b1_i = 0;
    for j = 1:i
        wk_return_b1_i = wk_return_b1_i + wk_return_b1(j);
    end
        sor_index(i) = sum(wk_return_b1_p(1:i)-rt_test(1:i)) / i / (dsr_index(i)^0.5);
        index_ac_wk_rt(i) = sum(wk_return_b1_p(1:i));
end


sor_yearly = sor * (290/5.5)^0.5;
sor_ew_yearly = sor_ew * (290/5.5)^0.5;
sor_index_yearly = sor_index * (290/5.5)^0.5;

%开始画图
%{
time = [1:145];
figure(3);
set(gcf, 'unit', 'centimeters', 'position', [10, 10, 16, 10])
marksize = 5;
linewidth_gca = 1.6;
fontsize_gca = 24;
fontsize_label = 26;
fontsize_legend = 20;
plot(time, sor_a_od(1:145), '-', 'linewidth', 1.6, 'markersize', marksize, 'color', [0 0.80784 0.81961]);
hold on
plot(time, sor(1:145), '--', 'linewidth', 1.6, 'markersize', marksize, 'color', '[0, 0, 0.9333]');
hold on
plot(time, sor_ew(1:145), '-.', 'linewidth', 1.6, 'markersize', marksize, 'color', 'r');
hold on
plot(time, sor_index(1:145), ':', 'linewidth', 1.6, 'markersize', marksize, 'color', 'm');

xlim([1 145])
ylim([0 1])
xticks([13 29 45 61 77 93 109 125 141])
xticklabels({'1995', '', '', '1996', '', '' ,'1997', '', ''})
h = legend('NO$_{\rm \bf a}$', 'NO$_{\rm \bf p}$', 'EW', 'MI');
set(h,'fontsize',fontsize_legend, 'interpreter', 'latex');
set(gca, 'linewidth', linewidth_gca, 'fontsize', fontsize_gca);
%set(gca, 'GridLineStyle', '--');
ylabel('SoR$_{\rm \bf p}$vs SoR$_{\rm \bf a}$', 'interpreter', 'latex', 'fontsize', fontsize_label)   ;
title('FTSE 100'); % 标题
hold on

hfig3 = figure(3);
figWidth = 16;  % 设置图片宽度
figHeight = 10;  % 设置图片高度
set(hfig3,'PaperUnits','centimeters'); % 图片尺寸所用单位
set(hfig3,'PaperPosition',[0 0 figWidth figHeight]);
%}


%各个sr的最终值比较
sor_yearly(end)
sor_ew_yearly(end)
sor_index_yearly(end)