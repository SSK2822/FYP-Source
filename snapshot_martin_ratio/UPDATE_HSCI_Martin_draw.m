clc,clear
format long;

filedir = './UPDATE_HSCI_Martin/';
filelist = dir([filedir, '*.txt']);
d1 = load('wk_price_HSCI', ',')';
b1 = load('b1_HSCI.txt', ','); %index
rf_92_97=load('rf_00_17.txt',','); %risk-free


[wk_return_d1, ~] = price2ret(d1', [], 'Periodic'); %simple return
wk_return_d1 = wk_return_d1';
[wk_return_b1, ~] = price2ret(b1, [], 'Periodic'); %simple return

[M, N] = size(wk_return_d1);
rf_test = rf_92_97(N/2+1:N);

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

for i = (N/2+1):N 
    wk_return_d1_train = wk_return_d1(:, 1:i-1); %training
    wk_return_d1_test = wk_return_d1(:, i); 
    rf = rf_92_97(i); 
    miu = mean(wk_return_d1_train, 2);
    filename = filelist(i-N/2);
    xt = load(['./UPDATE_HSCI_Martin/',filename.name]);
    xt_all(:,i-N/2) = xt;
    
    %martin_ratio_p
    My_wk_rt_temp = xt' * wk_return_d1_test - rf;
    My_wk_rt(i-N/2) = My_wk_rt_temp;
    martin_drawdown0 = Martin_Var_p(i-N/2, wk_return_d1, xt_all);
    ratio_my = sum(My_wk_rt) / (i-N/2) / (martin_drawdown0)^0.5;
    ratio_p_all(i-N/2) = ratio_my;

    %cumulative return
    ac_wk_rt_temp = ac_wk_rt_temp + xt' * wk_return_d1_test;
    ac_wk_rt(i-N/2) = ac_wk_rt_temp;

    %martin_ratio_a
    My_wk_rt_a_temp = xt' * mean(wk_return_d1_train, 2) - rf;
    ratio_a_my = My_wk_rt_a_temp / (Martin_Var_a(i, wk_return_d1, xt))^0.5;
    ratio_a_all(i-N/2) = ratio_a_my;
    
    %equally weighted
    ew_wk_rt_temp = x_ew' * wk_return_d1_test - rf;
    ew_wk_rt(i-N/2) = ew_wk_rt_temp;
    ratio_ew = sum(ew_wk_rt) / (i-N/2) / (Martin_Var_p(i-N/2, wk_return_d1, x_ew_all))^0.5;
    esr_ew(i-N/2) = ratio_ew;
    
    ew_ac_wk_rt_temp = ew_ac_wk_rt_temp + x_ew' * wk_return_d1_test;
    ew_ac_wk_rt(i-N/2) = ew_ac_wk_rt_temp;
    
end

%index怎么求
%index
esr_index = zeros(1, N/2);%初始化
index_wk_rt = 0;
for i = N/2+1:N
    index_wk_rt = index_wk_rt + wk_return_b1(:, i) - rf_92_97(:, i);
    %[~, evar_index, a] = fmincon(@(rho) EVaR_p_index(i, wk_return_b1, rho, theta), ...
    %10, [], [], [], [], 0, [], []);
    %disp(['index_p_cal:' num2str(a)]);
    evar_index = Martin_Var_p_index(i-N/2, wk_return_b1);
    esr_index(:, i-N/2) = index_wk_rt / (i-N/2) / (evar_index^0.5);
end

ratio_p_all_yearly = ratio_p_all * (N/18)^0.5;
ratio_a_all_yearly = ratio_a_all * (N/18)^0.5;
ac_wk_rt_yearly = ac_wk_rt * (N/18)^0.5;
esr_ew_yearly = esr_ew * (N/18)^0.5;
esr_index_yearly = esr_index * (N/18)^0.5;


%各个esr的最终值比较
ratio_p_all_yearly(end)
ratio_a_all_yearly(end)
ac_wk_rt_yearly(end)
esr_ew_yearly(end)
esr_index_yearly(end)