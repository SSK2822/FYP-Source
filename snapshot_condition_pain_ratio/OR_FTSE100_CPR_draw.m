clc,clear
format long;

filedir = './OR_FTSE100_CPR/';
filelist = dir([filedir, '*.txt']);
rhodir = './OR_FTSE100_CPR_rho/';
rholist = dir([rhodir, '*.txt']);
d1 = load('ORalldata_FTSE.txt', ',');
b1 = load('ORb1_FTSE.txt', ','); %index
rf_92_97=load('rf_92_97.txt',','); %risk-free


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
    i
    wk_return_d1_train = wk_return_d1(:, 1:i-1); %training
    wk_return_d1_test = wk_return_d1(:, i); 
    rf = rf_92_97(i); 
    miu = mean(wk_return_d1_train, 2);
    filename = filelist(i-N/2);
    xt = load(['./OR_FTSE100_CPR/',filename.name]);
    xt_all(:,i-N/2) = xt;
    
    %cpr_ratio_p
    My_wk_rt_temp = xt' * wk_return_d1_test - rf;
    My_wk_rt(i-N/2) = My_wk_rt_temp;
    cdd_drawdown0 = CDD_Var_p_cal(i-N/2, wk_return_d1, xt_all, theta);
    ratio_my = sum(My_wk_rt) / (i-N/2) / (cdd_drawdown0);
    ratio_p_all(i-N/2) = ratio_my;

    %cumulative return
    ac_wk_rt_temp = ac_wk_rt_temp + xt' * wk_return_d1_test;
    ac_wk_rt(i-N/2) = ac_wk_rt_temp;

    %cpr_ratio_a
    %{
    My_wk_rt_a_temp = xt' * mean(wk_return_d1_train, 2) - rf;
    ratio_a_my = My_wk_rt_a_temp / CDD_Var_a_cal(i, wk_return_d1, xt, theta);
    ratio_a_all(i-N/2) = ratio_a_my;
    %}

    %equally weighted
    ew_wk_rt_temp = x_ew' * wk_return_d1_test - rf;
    ew_wk_rt(i-N/2) = ew_wk_rt_temp;
    ratio_ew = sum(ew_wk_rt) / (i-N/2) / CDD_Var_p_cal(i-N/2, wk_return_d1, x_ew_all, theta);
    esr_ew(i-N/2) = ratio_ew;
    
    ew_ac_wk_rt_temp = ew_ac_wk_rt_temp + x_ew' * wk_return_d1_test;
    ew_ac_wk_rt(i-N/2) = ew_ac_wk_rt_temp;
    
end

%index怎么求 不需要x
%index
esr_index = zeros(1, N/2);%初始化
index_wk_rt = 0;
for i = N/2+1:N
    index_wk_rt = index_wk_rt + wk_return_b1(:, i) - rf_92_97(:, i);
    %[~, evar_index, a] = fmincon(@(rho) EVaR_p_index(i, wk_return_b1, rho, theta), ...
    %10, [], [], [], [], 0, [], []);
    %disp(['index_p_cal:' num2str(a)]);
    evar_index = CDD_Var_p_index_cal(i-N/2, wk_return_b1, theta);
    esr_index(:, i-N/2) = index_wk_rt / (i-N/2) / evar_index;
end

ratio_p_all_yearly = ratio_p_all * (N/5.5)^0.5;
ratio_a_all_yearly = ratio_a_all * (N/5.5)^0.5;
ac_wk_rt_yearly = ac_wk_rt * (N/5.5)^0.5;
esr_ew_yearly = esr_ew * (N/5.5)^0.5;
esr_index_yearly = esr_index * (N/5.5)^0.5;


%各个esr的最终值比较
ratio_p_all_yearly(end)
ratio_a_all_yearly(end)
ac_wk_rt_yearly(end)
esr_ew_yearly(end)
esr_index_yearly(end)