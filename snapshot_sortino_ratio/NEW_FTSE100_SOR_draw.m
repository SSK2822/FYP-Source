clc,clear
format long;

filedir = './BCST_FTSE100_SOR/';
filelist = dir([filedir, '*.txt']);
load('FTSE100.mat'); %stock price matrix(939 weeks)
wk_return_d1 = Assets_Returns(1:716, :)';
wk_return_b1 = Index_Returns(1:716, :)';
[M, N] = size(wk_return_d1);

rf_00_17=load('rf_02_16.txt',','); %risk free
rt_test = rf_00_17(N/2+1:N);%后半期risk-free

My_wk_rt_temp = 0; %初始化用于后面ratio分子的计算
ac_wk_rt_temp = 0;
ew_wk_rt_temp = 0;
ew_ac_wk_rt_temp = 0;
xt = zeros(M, 1);%初始化输出列向量xt
xt_all = zeros(M, N/2);
x_ew = 1 / M * ones(M, 1);
sor = zeros(1, N/2);%初始化行向量%目标输出测试集的所有sr
sor_ew = zeros(1, N/2);
My_wk_rt = zeros(1, N/2);
ac_wk_rt = zeros(1, N/2);
ew_wk_rt = zeros(1, N/2);
ew_ac_wk_rt = zeros(1, N/2);
sor_a_od = zeros(1, N/2);

for i = (N/2+1):N %后半期的第i周输出：xt，sor_p,绘图：
    wk_return_d1_train = wk_return_d1(:, 1:i-1); %用以训练的return
    wk_return_d1_test = wk_return_d1(:, i); %第i时期的return
    rf = rf_00_17(i); %第i时期的risk free
    miu = mean(wk_return_d1_train, 2); %每个时期的miu %每期更新
    filename = filelist(i-N/2);
    xt = load(['./BCST_FTSE100_SOR/',filename.name]);
    xt_all(:, i-N/2) = xt;
    
    %计算sor值  %可以顺便只求return
    My_wk_rt_temp = xt' * wk_return_d1_test - rf;
    My_wk_rt(i-N/2) = My_wk_rt_temp;
    ac_wk_rt_temp = ac_wk_rt_temp + xt' * wk_return_d1_test;
    ac_wk_rt(i-N/2) = ac_wk_rt_temp;
    dsr0 = DSR_p(i, wk_return_d1, xt_all, rf_00_17);
    ratio_my = sum(My_wk_rt) / (i-N/2) / (dsr0 ^ 0.5);
    sor(i-N/2) = ratio_my;

    %计算先验sor值
    My_wk_rt_a_temp = xt' * mean(wk_return_d1_train, 2) - rf;
    ratio_a_my = My_wk_rt_a_temp / (DSR_a(i, wk_return_d1, xt, rf_00_17))^0.5;
    sor_a_od(i-N/2) = ratio_a_my;

    %equally weighted sr_p
    ew_wk_rt_temp = x_ew' * wk_return_d1_test - rf;
    ew_wk_rt(i-N/2) = ew_wk_rt_temp;
    ew_ac_wk_rt_temp = ew_ac_wk_rt_temp + x_ew' * wk_return_d1_test;
    ew_ac_wk_rt(i-N/2) = ew_ac_wk_rt_temp;
    dsr0 = DSR_p_or(i-N/2, wk_return_d1, x_ew, rf_00_17);
    ratio_ew = sum(ew_wk_rt) / (i-N/2) / (dsr0 ^ 0.5);
    sor_ew(i-N/2) = ratio_ew;
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


sor_yearly = sor * (N/13.83)^0.5;
sor_ew_yearly = sor_ew * (N/13.83)^0.5;
sor_index_yearly = sor_index * (N/13.83)^0.5;


%各个sr的最终值比较
sor_yearly(end)
sor_ew_yearly(end)
sor_index_yearly(end)