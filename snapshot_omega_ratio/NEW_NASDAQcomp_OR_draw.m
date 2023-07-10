% clc,clear

filedir = './BCST_NASDAQcomp_OR/';
filelist = dir([filedir, '*.txt']);
load('NASDAQcomp.mat'); %stock price matrix(939 weeks)
wk_return_d1 = Assets_Returns(1:684, :)';
wk_return_b1 = Index_Returns(1:684, :)';
rf_90_16=load('rf_03_16.txt',','); %risk free
[M, N] = size(wk_return_d1);
rf_test = rf_90_16(N/2+1:N);


My_wk_rt_temp = 0; 
ac_wk_rt_temp = 0;
ew_wk_rt_temp = 0;
ew_ac_wk_rt_temp = 0;
xt = zeros(M, 1);
xt_all = zeros(M, N/2);
x_ew = 1 / M * ones(M, 1);
x_ew_all = 1 / M * ones(M, N/2);
or = zeros(1, N/2);
or_ew = zeros(1, N/2);
My_wk_rt = zeros(1, N/2);
ac_wk_rt = zeros(1, N/2);
ew_wk_rt = zeros(1, N/2);
ew_ac_wk_rt = zeros(1, N/2);
or_a = zeros(1, N/2);

%solutions

for i = (N/2+1):N %后半期的第i周输出：xt，sor_p：
    wk_return_d1_train = wk_return_d1(:, 1:i-1); %用以训练的return
    wk_return_d1_test = wk_return_d1(:, i); %第i时期的return
    rf = rf_90_16(i); %第i时期的risk free
    miu = mean(wk_return_d1_train, 2); %每个时期的miu %每期更新
    filename = filelist(i-N/2);
    xt = load(['./BCST_NASDAQcomp_OR/',filename.name]);
    xt_all(:, i-N/2) = xt;
    
    %or_p
    My_wk_rt_temp = xt' * wk_return_d1_test - rf;
    My_wk_rt(i-N/2) = My_wk_rt_temp;
    ac_wk_rt_temp = ac_wk_rt_temp + xt' * wk_return_d1_test;
    ac_wk_rt(i-N/2) = ac_wk_rt_temp;
    lpm0 = LPM_p(i, wk_return_d1, xt_all, rf_90_16);
    ratio_my = sum(My_wk_rt) / (i-N/2) / (lpm0);
    or(i-N/2) = ratio_my;

    %or_a
    My_wk_rt_a_temp = xt' * mean(wk_return_d1_train, 2) - rf;
    ratio_a_my = My_wk_rt_a_temp / (LPM_a(i, wk_return_d1, xt, rf_90_16));
    or_a(i-N/2) = ratio_a_my;

    %equally weighted sr_p
    ew_wk_rt_temp = x_ew' * wk_return_d1_test - rf;
    ew_wk_rt(i-N/2) = ew_wk_rt_temp;
    ew_ac_wk_rt_temp = ew_ac_wk_rt_temp + x_ew' * wk_return_d1_test;
    ew_ac_wk_rt(i-N/2) = ew_ac_wk_rt_temp;
    lpm0 = LPM_p(i, wk_return_d1, x_ew_all, rf_90_16);
    ratio_ew = sum(ew_wk_rt) / (i-N/2) / (lpm0);
    or_ew(i-N/2) = ratio_ew;
    
    
    
end


%任务1:sor的比较
%index的sor
wk_return_b1_p = wk_return_b1((N/2+1):N);
lpm_index = zeros(1, N/2); %初始化
lpm = zeros(1, N/2);%初始化
lpm_temp = zeros(1, N/2);
for i = 1:(N/2)
    for j = 1:i
        if (wk_return_b1_p(j) < rf_test(j))
            lpm(j) = (rf_test(j) - wk_return_b1_p(j));
            lpm_temp(i) = sum(lpm(1:i));
        end
        lpm_index(i) = lpm_temp(i) / i;
    end
end
or_index = zeros(1, N/2);
index_ac_wk_rt = zeros(1, N/2);
for i = 1:(N/2)
    or_index(i) = sum(wk_return_b1_p(1:i)-rf_test(1:i)) / i / (lpm_index(i));
    index_ac_wk_rt(i) = sum(wk_return_b1_p(1:i));
end

or_yearly = or * (N/13.25)^0.5;
or_ew_yearly = or_ew * (N/13.25)^0.5;
or_index_yearly = or_index * (N/13.25)^0.5;


%各个sor的最终值比较
or_yearly(end)
or_ew_yearly(end)
or_index_yearly(end)