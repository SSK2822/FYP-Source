clc,clear
format long;

file = './DowJones_SR55/';
load('DowJones.mat'); %stock price matrix(939 weeks)
wk_return_d1 = Assets_Returns(1:1362, :)';
wk_return_b1 = Index_Returns(1:1362, :)';
rf_90_16=load('rf_90_16.txt',','); %risk free

%variable initial
[M, N] = size(wk_return_d1);
rt_test = rf_90_16(N/2+1:N);
My_wk_rt_temp = 0; 
ac_wk_rt_temp = 0;
ew_wk_rt_temp = 0;
ew_ac_wk_rt_temp = 0;
xt = zeros(M, 1);
x_ew = 1 / M * ones(M, 1);
sr = zeros(1, N/2);
sr_ew = zeros(1, N/2);
My_wk_rt = zeros(1, N/2);
ac_wk_rt = zeros(1, N/2);
ew_wk_rt = zeros(1, N/2);
ew_ac_wk_rt = zeros(1, N/2);
sr_a = zeros(1, N/2);

for i = (N/2+1):N/2+1 
    wk_return_d1_train = wk_return_d1(:, 1:i-1); %training data
    wk_return_d1_test = wk_return_d1(:, i); %test data
    rf = rf_90_16(i); 
    miu = mean(wk_return_d1_train, 2); 
    if (miu < rf) 
        fprintf('do not trade in week %d\n', i)
        xt = zeros(M, 1); 
    else 
        xt = NNN_IG(i, wk_return_d1, rf_90_16, miu); %calculate decision xt using NO
    end 


end
