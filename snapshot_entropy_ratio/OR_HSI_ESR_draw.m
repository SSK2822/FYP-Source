clc,clear
format long;

filedir = './OR_HSI_ESR/';
filelist = dir([filedir, '*.txt']);
rhodir = './OR_HSI_ESR_rho/';
rholist = dir([rhodir, '*.txt']);
d1 = load('ORalldata_HS.txt', ',');
b1 = load('ORb1_HS.txt', ','); %index
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

for i = (N/2+1):N 
    wk_return_d1_train = wk_return_d1(:, 1:i-1); %training
    wk_return_d1_test = wk_return_d1(:, i); 
    rf = rf_92_97(i); 
    miu = mean(wk_return_d1_train, 2);
    filename = filelist(i-N/2);
    xt = load(['./OR_HSI_ESR/',filename.name]);
    xt_all(:,i-N/2) = xt;
    
    %ESR_p
    My_wk_rt_temp = xt' * wk_return_d1_test - rf;
    My_wk_rt(i-N/2) = My_wk_rt_temp;
    evar0 = EVaR_p_cal(i-N/2, wk_return_d1, xt_all, theta);
    % evar0 = EVaR_p(i-N/2, wk_return_d1, xt_all, theta);
    ratio_my = sum(My_wk_rt) / (i-N/2) / evar0;
    esr(i-N/2) = ratio_my;

    %cumulative return
    ac_wk_rt_temp = ac_wk_rt_temp + xt' * wk_return_d1_test;
    ac_wk_rt(i-N/2) = ac_wk_rt_temp;

    %esr_a
    if (i<N)
        filename_rho = rholist(i-N/2+1);
        rhot = load(['./OR_HSI_ESR_rho/',filename_rho.name]);
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


esr_yearly = esr * (N/5.5)^0.5;
esr_ew_yearly = esr_ew * (N/5.5)^0.5;


%开始画图

%各个esr的最终值比较
esr_yearly(end)
esr_ew_yearly(end)